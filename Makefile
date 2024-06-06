N = 26875598329180043462259934084516104814771548267173013831211000221401001686930205753736558686595756830038589223324899743441482269956074202215871484774010219580671726150483606643376601460289991046182995164525470357665094325803403228832463249478975783872562267277642071795087683984026446047408228788041596070593961009836625216432727028013369284809736590028172428593065991282443171883338828042258150663744334587327440242045135825821735601504456562981207961309937979063563597818822213630512868672556481966229905707148456095202031953667905249307250009435442291064206489495881279709522595506261753205200263886282673826508123

all: task

task: venv/bin/factorgpg gnupg.txt
	venv/bin/factorgpg gnupg.txt $(N)

venv/bin/factorgpg: venv/
	venv/bin/pip install -e .

venv/:
	python3 -m venv venv
	venv/bin/pip install --upgrade pip

clean:
	rm -rf build/
	rm -rf src/factorgpg/__pycache__
	rm -rf src/*.egg-info
	rm -rf venv/

.PHONY: all task clean
