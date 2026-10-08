# Pull GHS from PubChem
PubChem API calls are performed to receive GHS information by CID.

Custom run parameters are available at the top of each `"__main__"` call.

Run workflow
Install [uv](https://docs.astral.sh/uv/getting-started/installation/) for Python and Python dependency management.
```bash
		> cd python
		> uv sync
		> .venv\Scripts\activate
		> python "pulling ghs.py"
```