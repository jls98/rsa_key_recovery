N = 89269618190559985582507957643394420725636959471981573681053130515697627424575059088957619143562088220378551977725401056268570771814897717720008606718044585031642960387412397498609238190269222127555359410964757952284858157638473467630772405009202164869621956024705025971379492587060685995810780860921075070543779545079395366868765500822227304937206410029961875131892405411196221651587064806065493529825975204393347433929023221308899059489383911262814195269294290705256951464831506631555680723662332725245972309680900598977104331707977329766180202762607356658610721639420183615501550004397271275850116061925902479683595

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