import pathlib
import logging
import csv
from Bio import SeqIO, SeqRecord, Seq
from Bio.SeqFeature import SimpleLocation, CompoundLocation, SeqFeature
from pathlib import Path
from polars import DataFrame
import polars as pl

def seqfeature_hash(seqfeature: SeqFeature):
    """Hash for a SeqFeature objects"""
    return hash(str(seqfeature))

def assign_features_ssr(ssrdf: DataFrame, seqobj: Seq, circular: bool):
    if ssrdf.is_empty():
        ssrdf = ssrdf.with_columns(pl.lit(None).alias("annotations").cast(pl.List(pl.String)),
            pl.lit(None).alias("annotationobjects").cast(pl.List(pl.Int64)))
        return ssrdf

    # Test to see if dataframe has a position and length column
    features = sequence_to_features_df(seqobj, circular)
    ssrdf = ssrdf.with_columns((pl.col("start") + pl.col("repeat_len")*pl.col("count")-1).alias("end").cast(pl.Int32))

    # Annotation on Left edge
    left_edge = ssrdf.join_where(features, (pl.col("left_bound") <= pl.col("start"))
                                         & (pl.col("start") <= pl.col("right_bound"))
                                         & (pl.col("end") > pl.col("right_bound")))
    left_edge = left_edge.with_columns(pl.lit("edge").alias("featureclass"))

    # Annotation on right edge
    right_edge = ssrdf.join_where(features, (pl.col("left_bound") <= pl.col("end"))
                                         & (pl.col("end") <= pl.col("right_bound"))
                                         & (pl.col("start") < pl.col("left_bound")))
    right_edge = right_edge.with_columns(pl.lit("edge").alias("featureclass"))

    # Annotation inside
    inside = ssrdf.join_where(features, (pl.col("start") <= pl.col("left_bound")) & (pl.col("end") >= pl.col("right_bound")))
    inside = inside.with_columns(pl.lit("inside").alias("featureclass"))

    # Annotation wraps
    wraps = ssrdf.join_where(features, (pl.col("start") >= pl.col("left_bound")) & (pl.col("end") <= pl.col("right_bound")))
    wraps = wraps.with_columns(pl.lit("wraps").alias("featureclass"))

    anno = pl.concat([left_edge, right_edge, inside, wraps]).unique()
    anno = anno.group_by(["repeat", "repeat_len", "start", "count", "mutation_rate"]
                        ).agg(pl.col('annotations'), pl.col('annotationobjects')
    )

    intergenic = ssrdf.join(anno.select(["repeat", "repeat_len", "start", "count"]), how="anti", on=["repeat", "repeat_len", "start", "count"]
    ).with_columns(pl.lit([]).alias("annotations")).with_columns(pl.lit([]).alias("annotationobjects")).select(["repeat", "repeat_len", "start", "count", "mutation_rate", "annotations", "annotationobjects"])


    ssrdf = pl.concat([anno, intergenic])

    return ssrdf


def assign_features_rmd(rmd_or_srs_df: DataFrame, seqobj: Seq, circular: bool):

    if rmd_or_srs_df.is_empty():
        rmd_or_srs_df = rmd_or_srs_df.with_columns(pl.lit(None).alias("annotations").cast(pl.List(pl.String)),
            pl.lit(None).alias("annotationobjects").cast(pl.List(pl.Int64)))
        return rmd_or_srs_df

    # Test to see if a dataframe has a position left, position right, and length column
    features = sequence_to_features_df(seqobj, circular)
    pl.Config.set_tbl_rows(100)
    sequence_length = len(seqobj)
    df = (rmd_or_srs_df
            .with_columns(pl.when(
            pl.col("second_repeat")-pl.col("first_repeat")-pl.col("repeat_len") != pl.col("distance")
            ).then(pl.lit(True)
            ).otherwise(pl.lit(False)
            ).alias("wraps")))

    df =    (df.with_columns(pl.when(
            pl.col("wraps") == False
            ).then(pl.col("first_repeat"))
            .otherwise(pl.col("second_repeat")).alias("start")
            )
            .with_columns(
            pl.when(
                pl.col("wraps") == False
            ).then(pl.col("second_repeat") + pl.col("repeat_len"))
            .otherwise(pl.col("first_repeat") + pl.col("repeat_len")).alias("end")
            ))

    # Annotation on Left edge
    left_edge = df.filter(pl.col("wraps") == False).join_where(features, ((pl.col("left_bound") <= pl.col("start")))
                                            & (pl.col("start") <= pl.col("right_bound"))
                                            & (pl.col("end") > pl.col("right_bound"))
                                            )
    left_edge_wraps = df.filter(pl.col("wraps") == True).join_where(features, (
                                            pl.col("left_bound") <= pl.col("start"))
                                            & (pl.col("right_bound") >= pl.col("start")))

    left_edge = pl.concat([left_edge, left_edge_wraps]).with_columns(pl.lit("edge").alias("featureclass"))

    # Annotation on right edge
    right_edge = df.filter(pl.col("wraps") == False).join_where(features, (pl.col("left_bound") <= pl.col("end"))
                                            & (pl.col("end") <= pl.col("right_bound"))
                                            & (pl.col("start") < pl.col("left_bound")))
    right_edge_wraps = df.filter(pl.col("wraps") == True).join_where(features, (
                                            pl.col("left_bound") <= pl.col("end"))
                                            & (pl.col("right_bound") >= pl.col("end")))
    right_edge = pl.concat([right_edge, right_edge_wraps]).with_columns(pl.lit("edge").alias("featureclass"))

    # Annotation inside
    inside = df.filter(pl.col("wraps") == False).join_where(features, (pl.col("start") <= pl.col("left_bound")) & (pl.col("end") >= pl.col("right_bound")))
    inside_wraps_a = df.filter(pl.col("wraps") == True).join_where(features, (pl.col("start") <= pl.col("left_bound")))
    inside_wraps_b = df.filter(pl.col("wraps") == True).join_where(features, (pl.col("end") >= pl.col("right_bound")))
    inside = pl.concat([inside, inside_wraps_a, inside_wraps_b]).with_columns(pl.lit("inside").alias("featureclass"))

    # Annotation wraps
    wraps = df.filter(pl.col("wraps") == False).join_where(features, (pl.col("start") >= pl.col("left_bound")) & (pl.col("end") <= pl.col("right_bound")))
    wraps = pl.concat([wraps]).with_columns(pl.lit("wraps").alias("featureclass"))

    anno = pl.concat([left_edge, right_edge, inside, wraps]).unique()
    anno = anno.group_by(["repeat", "repeat_len", "first_repeat", "second_repeat", "distance", "mutation_rate"]
                        ).agg(pl.col('annotations'), pl.col('annotationobjects'))

    intergenic = df.join(anno.select(["repeat", "repeat_len", "first_repeat", "second_repeat", "distance"]), how="anti", on=["repeat", "repeat_len", "first_repeat", "second_repeat", "distance"]
    ).with_columns(pl.lit([]).alias("annotations")).with_columns(pl.lit([]).alias("annotationobjects")).select(["repeat", "repeat_len", "first_repeat", "second_repeat", "distance", "mutation_rate", "annotations", "annotationobjects"])

    df = pl.concat([anno, intergenic])

    return df


FASTA_EXTS = [".fa", ".fasta"]
GBK_EXTS = [".gb", ".gbk", ".gbff"]

def determine_looparound(df):
    """Determines if the prediction wraps around the origin"""
    pass

def sequence_to_features_df(sequence, circular=True):
    """Takes genbank annotations and turns them into a polars dataframe"""
    features = sequence.features
    seqlen = len(sequence)

    def get_feature_bounds(feature_location):
        if isinstance(feature_location, SimpleLocation):
            return feature_location.start, feature_location.end
        elif isinstance(feature_location, CompoundLocation):
            if feature_location.parts[-1].end < feature_location.parts[0].start:
                end = len(sequence) + feature_location.parts[-1].end
            else:
                end = feature_location.parts[-1].end
            return feature_location.parts[0].start, end

    if not circular:
        # Look for compound wraparounds and break them into simple features
        newfeatures = []
        deletedfeatures = []
        for feature in features:
            if not isinstance(feature.location, CompoundLocation):
                continue

            # Check whether the compound feature is actually a wraparound
            wraparound_part_index = None
            last_part_start = None
            rightmost_part = None
            for i, part in enumerate(feature.location.parts):
                if rightmost_part != None and part.start < last_part_start:
                    wraparound_part_index = i
                if rightmost_part == None:
                    rightmost_part = i
                    last_part_start = part.start
            if wraparound_part_index is None:
                continue

            # If it is a wraparound, break it into two features
            leftsplit_locations = feature.location.parts[:wraparound_part_index]
            if len(leftsplit_locations) > 1:
                leftsplit_locations = CompoundLocation(leftsplit_locations)
            else:
                leftsplit_locations = leftsplit_locations[0]
            newleftfeature = SeqFeature(leftsplit_locations, feature.type, qualifiers=feature.qualifiers)


            rightsplit_locations = feature.location.parts[wraparound_part_index:]
            if len(rightsplit_locations) > 1:
                rightsplit_locations = CompoundLocation(rightsplit_locations)
            else:
                rightsplit_locations = rightsplit_locations[0]
            newrightfeature = SeqFeature(rightsplit_locations, feature.type, qualifiers=feature.qualifiers)

            newfeatures.append(newleftfeature)
            newfeatures.append(newrightfeature)
            deletedfeatures.append(feature)

        for feature in deletedfeatures:
            features.remove(feature)
        features.extend(newfeatures)

    df = pl.DataFrame([(feature.type,
                        get_feature_bounds(feature.location),
                        feature.qualifiers.get("label", ["unlabeled"])[0],
                        seqfeature_hash(feature)) for feature in features],
        schema=['type', 'loc', 'annotations', 'annotationobjects'],
        orient="row")

    # expand out loc
    df = df.with_columns(pl.col("loc").list.to_struct(fields=['left_bound', 'right_bound'])).unnest("loc")
    df = df.with_columns(
        pl.col("left_bound").cast(pl.Int32),
        pl.col("right_bound").cast(pl.Int32)
    )

    return df
