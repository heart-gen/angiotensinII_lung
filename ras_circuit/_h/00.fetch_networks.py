"""
Fetch and cache the decoupler prior networks used by 03.at1r_response_score.py.

PROGENy (pathway footprints) and CollecTRI (TF regulons) are downloaded once,
through OmniPath, and written to _m/networks/ with a PROVENANCE.txt. Every later
step reads ONLY the cached files, so a rerun cannot silently pick up a newer
network release. Delete _m/networks/ deliberately to refresh.

Compute nodes have internet access here, so this runs as the first command of
step_2.sh rather than as a separate login-node step.
"""
import argparse
import hashlib
import logging
from datetime import datetime, timezone
from pathlib import Path

import decoupler as dc
import omnipath


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def with_retries(fn, what, tries=4, wait=20):
    """Zenodo occasionally drops the TLS stream mid-download (seen 2026-09-22:
    SSLEOFError from a compute node, fine from the login node minutes later)."""
    import time
    for k in range(1, tries + 1):
        try:
            return fn()
        except Exception as e:  # network errors surface as several types
            logging.warning(f"{what}: attempt {k}/{tries} failed ({type(e).__name__}: {e})")
            if k == tries:
                raise
            time.sleep(wait * k)


def collectri_via_omnipath():
    """Fallback route: the same CollecTRI resource served by OmniPath, converted to
    decoupler's (source, target, weight) form -- weight -1 for inhibitory
    consensus, +1 otherwise, which is how decoupler itself signs the network."""
    ct = omnipath.interactions.CollecTRI.get(genesymbols=True, organism="human", loops=True)
    w = ct["consensus_inhibition"].astype(bool) & ~ct["consensus_stimulation"].astype(bool)
    net = ct.assign(source=ct["source_genesymbol"], target=ct["target_genesymbol"],
                    weight=(-1.0) * w + 1.0 * ~w)[["source", "target", "weight"]]
    return net.drop_duplicates(["source", "target"])


def main():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--outdir", type=Path, default=Path("./networks"))
    p.add_argument("--progeny-top", type=int, default=500)
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    args.outdir.mkdir(parents=True, exist_ok=True)

    targets = {
        f"progeny_human_top{args.progeny_top}.tsv":
            lambda: dc.op.progeny(organism="human", top=args.progeny_top),
        "collectri_human.tsv": lambda: dc.op.collectri(organism="human"),
    }
    lines = [f"fetched_utc\t{datetime.now(timezone.utc).isoformat()}",
             f"decoupler\t{dc.__version__}",
             f"omnipath\t{omnipath.__version__}"]
    for fname, fetch in targets.items():
        out = args.outdir / fname
        if out.exists():
            logging.info(f"{fname}: cached, not re-fetched")
        else:
            route = "decoupler.op"
            try:
                net = with_retries(fetch, fname)
            except Exception:
                if not fname.startswith("collectri"):
                    raise
                logging.warning("CollecTRI: decoupler route failed; using OmniPath directly")
                net = with_retries(collectri_via_omnipath, "CollecTRI via OmniPath")
                route = "omnipath.interactions.CollecTRI (fallback)"
            lines.append(f"{fname}\troute={route}")
            net.to_csv(out, sep="\t", index=False)
            logging.info(f"{fname}: {net.shape[0]} edges, "
                         f"{net['source'].nunique()} sources")
        lines.append(f"{fname}\tmd5={md5(out)}")
    prov = args.outdir / "PROVENANCE.txt"
    if not prov.exists():
        prov.write_text("\n".join(lines) + "\n")
    logging.info("networks ready in " + str(args.outdir))


if __name__ == "__main__":
    main()
