raise NotImplementedError()
import pathlib
import logging
import csv
from Bio import SeqIO, SeqRecord, Seq
from Bio.SeqFeature import SimpleLocation, CompoundLocation, SeqFeature
from pathlib import Path
from polars import DataFrame
import polars as pl


def assign_features_ssr(ssrdf: DataFrame, seqobj: Seq, circular: bool):
    # Test to see if dataframe has a position and length column
    features = sequence_to_features_df(seqobj, circular)
    ssrdf = ssrdf.with_columns((pl.col("start") + pl.col("repeat_len")*pl.col("count")-1).alias("end"))

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
                        ).agg(pl.col('name'), pl.col('object')
    )

    intergenic = ssrdf.join(anno.select(["repeat", "repeat_len", "start", "count"]), how="anti", on=["repeat", "repeat_len", "start", "count"]
    ).with_columns(pl.lit([]).alias("name")).with_columns(pl.lit([]).alias("object")).select(["repeat", "repeat_len", "start", "count", "mutation_rate", "name", "object"])


    ssrdf = pl.concat([anno, intergenic])

    return ssrdf


def assign_features_rmd(rmd_or_ssr_df: DataFrame, seqobj: Seq, circular: bool):
    # Test to see if a dataframe has a position left, position right, and length column
    features = sequence_to_features_df(seqobj, circular)
    sequence_length = len(seqobj)
    df = (rmd_or_ssr_df
            .with_columns(pl.when(
            pl.col("distance") > sequence_length
            ).then(pl.lit(True)
            ).otherwise(pl.lit(False)
            ).alias("wraps"))

            .with_columns(pl.when(
            pl.col("wraps") == False
            ).then(pl.col("position_left"))
            .otherwise(pl.col("position_right")).alias("start")
            )
            .with_columns(
            pl.when(
                pl.col("wraps") == False
            ).then(pl.col("position_right") + pl.col("repeat_len"))
            .otherwise(pl.col("position_left") + pl.col("repeat_len")).alias("end")
            ))

    # ----- Non-wrapping examples

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
    wraps_wraps_a = df.filter(pl.col("wraps") == True).join_where(features, (pl.col("start") >= pl.col("left_bound")))
    wraps_wraps_b = df.filter(pl.col("wraps") == True).join_where(features, (pl.col("end") <= pl.col("right_bound")))
    wraps = pl.concat([wraps, wraps_wraps_a, wraps_wraps_b]).with_columns(pl.lit("wraps").alias("featureclass"))

    anno = pl.concat([left_edge, right_edge, inside, wraps]).unique()
    anno = anno.group_by(["repeat", "repeat_len", "position_left", "position_right", "distance", "mutation_rate"]
                        ).agg(pl.col('name'), pl.col('object'))

    intergenic = df.join(anno.select(["repeat", "repeat_len", "position_left", "position_right", "distance"]), how="anti", on=["repeat", "repeat_len", "position_left", "position_right", "distance"]
    ).with_columns(pl.lit([]).alias("name")).with_columns(pl.lit([]).alias("object")).select(["repeat", "repeat_len", "position_left", "position_right", "distance", "mutation_rate", "name", "object"])

    df = pl.concat([anno, intergenic])

    return df


def assign_features():
    pass


def position_assigner(seqobj):
    def inner(positions: tuple[Int32], circular: bool):
        pass
    return assign_feature

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
                        feature.qualifiers.get("label", "")[0],
                        feature) for feature in features],
        schema=['type', 'loc', 'name', 'object'],
        orient="row")

    # expand out loc
    df = df.with_columns(pl.col("loc").list.to_struct(fields=['left_bound', 'right_bound'])).unnest("loc")

    return df


def parse_file(filepath: pathlib.Path, use_filename: bool = True) -> list:
    """
    parses the inputted files and returns a list of sequences found in each file

    input:
        filepath: path to the multifasta, genbank, or csv file containing the sequences to be scanned

    returns:
        list containing all the sequences to be scanned
    """

    path_as_string = str(filepath)
    if not filepath.exists():
        raise OSError("File {} does not exist.".format(path_as_string))
    elif filepath.suffix in FASTA_EXTS:
        sequences = SeqIO.parse(path_as_string, "fasta")
    elif filepath.suffix in GBK_EXTS:
        sequences = SeqIO.parse(path_as_string, "genbank")
    else:
        raise ValueError(
            f"File {filepath} is not a known file format. Must be one of {FASTA_EXTS + GBK_EXTS + [".csv"]}."
        )

    # If the genbank file doesnt have a name, add it
    test_sequences = []
    for i, seq in enumerate(sequences):
        try:
            if use_filename:
                filename = Path(filepath).stem
                if not seq.name:
                    seq.name = f"{filename}"
                if not seq.description or seq.description == '':
                    seq.description = f"{filename}"
            test_sequences.append(seq)
        except:
            pass

    return test_sequences

if __name__ == "__main__":
    rmd_df = pl.read_csv("rmd_df")
    ssr_df = pl.read_csv("ssr_df")
    srs_df = pl.read_csv("srs_df")

    sequence = parse_file(Path("l6-10_plasmid_bba.gb"))[0]
    assign_features_ssr(ssr_df, sequence, circular=True)
