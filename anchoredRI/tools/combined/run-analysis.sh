set -eu
python3 audit.py
R --version | head -1
python3 score.py
python3 verify-scores.py
python3 audit-routes.py
python3 report.py
python3 write-report.py
cat results/REPORT.md
