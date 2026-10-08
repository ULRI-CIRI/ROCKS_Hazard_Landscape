import os
import time
import shutil
import socket
import pandas as pd
import requests

# ------------------------------------------------------------------
# FORCE REQUEST TIMEOUTS
# ------------------------------------------------------------------
# This monkey-patches requests so that if pubchem.py uses requests.get()
# without a timeout, it still gets a timeout.
#
# connect timeout = 10 sec
# read timeout    = 60 sec

_original_request = requests.sessions.Session.request

def request_with_timeout(self, method, url, **kwargs):
    kwargs.setdefault("timeout", (10, 60))
    return _original_request(self, method, url, **kwargs)

requests.sessions.Session.request = request_with_timeout

from pubchem import get_ghs_data

# ------------------------------------------------------------------
# INPUTS
# ------------------------------------------------------------------

cid_file = r"c:\projects\toxicity_landscape\data\raw\cids.txt"

outfile = (
    r"c:\projects\toxicity_landscape\data\derived\ghs_results.csv"
)

save_every = 500
progress_every = 100

max_attempts = 5
initial_retry_sleep = 10
request_pause = 0.25

max_consecutive_transient_errors = 25

pubchem_host = "pubchem.ncbi.nlm.nih.gov"

output_columns = [
    "inchikey",
    "cid",
    "error",
    "compound_name",
    "signal_word",
    "hazard_statements",
    "attempts"
]

# ------------------------------------------------------------------
# HELPER FUNCTIONS
# ------------------------------------------------------------------

def is_transient_error(x):
    """
    Return True for errors that should be retried on a future run.
    These usually reflect network, DNS, timeout, or connection issues,
    not true absence of GHS data.
    """

    if pd.isna(x):
        return False

    msg = str(x)

    transient_patterns = [
        "NameResolutionError",
        "Failed to resolve",
        "getaddrinfo failed",
        "Max retries exceeded",
        "HTTPSConnectionPool",
        "ConnectionError",
        "ConnectTimeout",
        "ReadTimeout",
        "Timeout",
        "timed out",
        "Temporary failure",
        "Connection aborted",
        "RemoteDisconnected",
        "ProtocolError",
        "ProxyError",
        "SSLError",
        "503",
        "502",
        "504"
    ]

    return any(pat in msg for pat in transient_patterns)


def check_dns(host):
    """
    Check whether the PubChem hostname resolves before starting.
    """

    try:
        socket.gethostbyname(host)
        return True
    except Exception:
        return False


def backup_existing_file(path):
    """
    Make a timestamped backup of the existing results file.
    """

    if not os.path.exists(path):
        return None

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    backup_path = path.replace(".csv", f"_backup_{timestamp}.csv")

    shutil.copy2(path, backup_path)

    print(
        f"Backed up existing output to:\n{backup_path}",
        flush=True
    )

    return backup_path


def clean_existing_results(path):
    """
    Remove retryable/transient network-error rows from the existing file.
    Keep completed rows and non-transient error rows.

    This prevents old DNS-failure rows from being treated as completed.
    """

    if not os.path.exists(path):
        return pd.DataFrame(columns=output_columns), set()

    try:
        existing = pd.read_csv(path)

    except Exception as e:
        print(
            f"Could not read existing output file: {e}",
            flush=True
        )
        return pd.DataFrame(columns=output_columns), set()

    # Make sure all expected columns exist
    for col in output_columns:
        if col not in existing.columns:
            existing[col] = None

    existing = existing[output_columns].copy()

    existing["cid"] = pd.to_numeric(
        existing["cid"],
        errors="coerce"
    )

    existing = existing.dropna(subset=["cid"]).copy()
    existing["cid"] = existing["cid"].astype(int)

    existing["retryable_error"] = existing["error"].apply(
        is_transient_error
    )

    n_total = len(existing)
    n_retryable = int(existing["retryable_error"].sum())

    cleaned = existing[
        ~existing["retryable_error"]
    ].copy()

    cleaned = cleaned.drop(columns=["retryable_error"])

    # If a CID appears multiple times, keep the last non-retryable row
    cleaned = cleaned.drop_duplicates(
        subset=["cid"],
        keep="last"
    )

    completed = set(cleaned["cid"].astype(int))

    # Rewrite outfile with cleaned rows only
    cleaned.to_csv(
        path,
        index=False
    )

    print(
        f"Existing rows before cleaning: {n_total:,}",
        flush=True
    )

    print(
        f"Removed retryable network/DNS error rows: {n_retryable:,}",
        flush=True
    )

    print(
        f"Completed CIDs retained: {len(completed):,}",
        flush=True
    )

    return cleaned, completed


def save_buffer(buffer, path):
    """
    Append buffered records to disk.
    """

    if not buffer:
        return

    pd.DataFrame(buffer).to_csv(
        path,
        mode="a",
        header=(
            not os.path.exists(path)
            or os.path.getsize(path) == 0
        ),
        index=False
    )


def get_ghs_with_retries(cid):
    """
    Try to retrieve GHS data for one CID.

    Returns:
        record dictionary
        transient_failure boolean
    """

    last_error = None

    for attempt in range(1, max_attempts + 1):

        try:
            ghs, err = get_ghs_data(int(cid))

            record = {
                "cid": int(cid),
                "error": err,
                "compound_name": getattr(
                    ghs,
                    "compound_name",
                    None
                ),
                "signal_word": getattr(
                    ghs,
                    "signal_word",
                    None
                ),
                "hazard_statements": (
                    "; ".join(ghs.hazard_statements)
                    if getattr(ghs, "hazard_statements", None)
                    else None
                ),
                "attempts": attempt
            }

            return record, False

        except Exception as e:

            last_error = str(e)
            transient = is_transient_error(last_error)

            print(
                (
                    f"CID {cid} failed on attempt "
                    f"{attempt}/{max_attempts}: {last_error}"
                ),
                flush=True
            )

            if attempt < max_attempts and transient:

                sleep_time = initial_retry_sleep * attempt

                print(
                    f"Retrying CID {cid} after {sleep_time} seconds...",
                    flush=True
                )

                time.sleep(sleep_time)

            else:

                record = {
                    "cid": int(cid),
                    "error": last_error,
                    "compound_name": None,
                    "signal_word": None,
                    "hazard_statements": None,
                    "attempts": attempt
                }

                return record, transient

    record = {
        "cid": int(cid),
        "error": last_error,
        "compound_name": None,
        "signal_word": None,
        "hazard_statements": None,
        "attempts": max_attempts
    }

    return record, is_transient_error(last_error)


# ------------------------------------------------------------------
# STARTUP CHECK
# ------------------------------------------------------------------

print(
    "Starting PubChem GHS pull...",
    flush=True
)

if not check_dns(pubchem_host):

    print(
        (
            f"DNS check failed for {pubchem_host}.\n"
            "PubChem hostname is not resolving right now. "
            "Stopping before writing more error rows."
        ),
        flush=True
    )

    raise SystemExit

print(
    f"DNS check passed for {pubchem_host}",
    flush=True
)

# ------------------------------------------------------------------
# LOAD CID LIST
# ------------------------------------------------------------------

cids = pd.read_csv(
    cid_file,
    sep="\t",
    names=["inchikey", "cid"]
)

cids["cid"] = pd.to_numeric(
    cids["cid"],
    errors="coerce"
)

cids = (
    cids
    .dropna(subset=["cid"])
    .drop_duplicates(subset=["cid"])
    .copy()
)

cids["cid"] = cids["cid"].astype(int)

print(
    f"Loaded {len(cids):,} unique CIDs",
    flush=True
)

# ------------------------------------------------------------------
# BACKUP AND CLEAN PRIOR OUTPUT
# ------------------------------------------------------------------

backup_existing_file(outfile)

_, completed = clean_existing_results(outfile)

remaining = cids[
    ~cids["cid"].isin(completed)
].copy()

print(
    f"{len(remaining):,} CIDs remaining",
    flush=True
)

if len(remaining) == 0:
    print("Nothing left to process.", flush=True)
    raise SystemExit

# ------------------------------------------------------------------
# PROCESS
# ------------------------------------------------------------------

buffer = []
consecutive_transient_errors = 0

for i, row in enumerate(
    remaining.itertuples(index=False),
    start=1
):

    inchikey = row.inchikey
    cid = int(row.cid)

    if i == 1 or i % progress_every == 0:

        print(
            (
                f"Processing CID {cid:,} | "
                f"{i:,}/{len(remaining):,} "
                f"remaining records processed"
            ),
            flush=True
        )

    try:

        record, transient_failure = get_ghs_with_retries(cid)

        record["inchikey"] = inchikey

        # Reorder columns
        record = {
            col: record.get(col, None)
            for col in output_columns
        }

        buffer.append(record)

        if transient_failure:
            consecutive_transient_errors += 1
        else:
            consecutive_transient_errors = 0

        if consecutive_transient_errors >= max_consecutive_transient_errors:

            print(
                (
                    "\nToo many consecutive transient network/DNS errors.\n"
                    f"Consecutive transient errors: "
                    f"{consecutive_transient_errors}\n"
                    "Saving current buffer and stopping to avoid filling "
                    "the CSV with network-failure rows."
                ),
                flush=True
            )

            save_buffer(buffer, outfile)
            buffer = []

            raise SystemExit

    except KeyboardInterrupt:

        print(
            "\nInterrupted by user. Saving current buffer...",
            flush=True
        )

        save_buffer(buffer, outfile)
        buffer = []

        raise

    # --------------------------------------------------------------
    # SAVE CHECKPOINT
    # --------------------------------------------------------------

    if len(buffer) >= save_every:

        save_buffer(buffer, outfile)

        print(
            f"Saved checkpoint ({len(buffer):,} records)",
            flush=True
        )

        buffer = []

    # Pause slightly to reduce request pressure
    time.sleep(request_pause)

# ------------------------------------------------------------------
# FINAL SAVE
# ------------------------------------------------------------------

if buffer:

    save_buffer(buffer, outfile)

    print(
        f"Saved final buffer ({len(buffer):,} records)",
        flush=True
    )

print(
    "Finished successfully",
    flush=True
)