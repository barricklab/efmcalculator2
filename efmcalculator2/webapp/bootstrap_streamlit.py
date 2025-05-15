import subprocess

from efmcalculator2 import run_webapp


def bootstrap_streamlit():
    subprocess.run(["streamlit", "run", __file__])


if __name__ == "__main__":
    run_webapp()
