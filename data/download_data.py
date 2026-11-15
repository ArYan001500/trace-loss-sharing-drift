"""Downloads the two public datasets used in the paper and places the files where the scripts expect them.

1. Traction-rig recordings, Zenodo record 20395811 (CC BY 4.0). Two ZIP archives are downloaded and checked
   against their published MD5 sums; only the four files used here are extracted to
   02_Data/raw/zenodo_20395811/extracted/.
2. Buck-converter thermal dataset, github.com/jm-sanz/Automatic-Thermal-Modeling (MIT license), pinned to commit
   d98fc094. Dataset.rar is downloaded and extracted, the four CSV files are checked against their MD5 sums, and
   make_sanz_1s.py builds 02_Data/working/sanz_1s/sanz_1s.npz (1 s means).

Extracting the RAR archive needs one of: bsdtar (the tar of Windows 10 and later, or libarchive-tools on Linux),
7z, or unrar.
python data/download_data.py"""
import hashlib, os, shutil, subprocess, sys, tempfile, urllib.request, zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RIG = os.path.join(ROOT, "02_Data", "raw", "zenodo_20395811", "extracted")
BOARD = os.path.join(ROOT, "02_Data", "working", "sanz_1s")

ZENODO = "https://zenodo.org/api/records/20395811/files/{}/content"
ZIPS = {
    "Mission profile RAW data.zip": ("17eea94180c904cf559c8d927baefe02",
                                     ["Mission profile/fr06final.csv", "Mission profile/fr06final_temperature.xlsx"]),
    "Parallel modules under different cooling circuits RAW data.zip": (
        "21e0e2cc6e345ed63a674982197913aa",
        ["Parallel modules under different cooling circuits MAG_1_1/L0_0/MAG1_1_L0_0.xlsx",
         "Parallel modules under different cooling circuits MAG_1_1/L0_50/MAG1_1_L0_50.xlsx"]),
}
RAR_URL = "https://github.com/jm-sanz/Automatic-Thermal-Modeling/raw/d98fc094da571536fc95a68d18e7a4418eb2bc48/Dataset.rar"
CSV_MD5 = {"Temp_Powers_Sat.csv": "a23f7b5af7e2581984e9d50100f6a930", "Temp_Powers_Lin.csv": "dfaf29a762e92ba0cbe5628726dcb5de",
           "Temp_Powers_Drv.csv": "04de7e70a1e9ac862526d5270f25dc63", "Temp_Powers_Ind.csv": "d035b190a293d9ed273c0dddb8faff78"}


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fetch(url, dest):
    print("downloading", url, flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "dataset-downloader"})
    with urllib.request.urlopen(req) as r, open(dest, "wb") as f:
        shutil.copyfileobj(r, f)


def rig(tmp):
    need = [n for _, names in ZIPS.values() for n in names]
    if all(os.path.exists(os.path.join(RIG, n)) for n in need):
        print("rig files present"); return
    for name, (checksum, members) in ZIPS.items():
        z = os.path.join(tmp, "rig.zip")
        fetch(ZENODO.format(urllib.request.quote(name)), z)
        if md5(z) != checksum:
            sys.exit(f"MD5 mismatch for {name}")
        with zipfile.ZipFile(z) as zf:
            for m in members:
                out = os.path.join(RIG, m)
                os.makedirs(os.path.dirname(out), exist_ok=True)
                with zf.open(m) as src, open(out, "wb") as dst:
                    shutil.copyfileobj(src, dst)
        os.remove(z)
    print("rig files extracted to", RIG)


def extract_rar(rar, dest):
    for cmd in (["bsdtar", "-xf", rar, "-C", dest], ["tar", "-xf", rar, "-C", dest],
                ["7z", "x", "-y", f"-o{dest}", rar], ["unrar", "x", "-o+", rar, dest + os.sep]):
        if shutil.which(cmd[0]) is None:
            continue
        if subprocess.run(cmd, capture_output=True).returncode == 0 and any(
                f.endswith(".csv") for _, _, fs in os.walk(dest) for f in fs):
            return
    sys.exit("Could not extract Dataset.rar: install bsdtar (libarchive-tools), 7z or unrar.")


def board(tmp):
    if os.path.exists(os.path.join(BOARD, "sanz_1s.npz")):
        print("board file present"); return
    rar = os.path.join(tmp, "Dataset.rar"); out = os.path.join(tmp, "board")
    os.makedirs(out)
    fetch(RAR_URL, rar)
    extract_rar(rar, out)
    found = {f: os.path.join(r, f) for r, _, fs in os.walk(out) for f in fs if f in CSV_MD5}
    for f, checksum in CSV_MD5.items():
        if f not in found or md5(found[f]) != checksum:
            sys.exit(f"missing or changed board file {f}")
    os.makedirs(BOARD, exist_ok=True)
    subprocess.run([sys.executable, os.path.join(ROOT, "data", "make_sanz_1s.py"), os.path.dirname(found["Temp_Powers_Sat.csv"]),
                    os.path.join(BOARD, "sanz_1s.npz")], check=True)
    print("board file written to", BOARD)


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmp:
        rig(tmp)
        board(tmp)
