date
source setup_main_injector.sh
python checkevent.py
curl http://textbelt.com/text -d number=2153008763 -d "message=CHECKEVENT.PY ENDED! MUST RESTART"
echo ""
date
