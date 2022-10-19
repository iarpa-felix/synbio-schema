Ubuntu 20. System has run similar projects in the past so can't be considers a fresh install.

- clone https://github.com/semantic-synbio/synbio-schema
- ensure your computer has python 3.9 as the default python on the system path ?
- ensure your system has the python poetry application installed and on the default path

```shell
python --version
```

> Python 3.9.5

```shell
poetry --version
```

> Poetry version 1.1.11

possibly delete poetry.lock ?

```shell
poetry install
```

> The Poetry configuration is invalid: Additional properties are not allowed ('group' was unexpected)

```
[[package]]
name = "linkml"
version = "1.3.8"
description = "Linked Open Data Modeling Language"
category = "dev"
optional = false
python-versions = ">=3.7.6,<4.0.0"

[[package]]
name = "linkml-runtime"
version = "1.3.4"
description = "Runtime environment for LinkML, the Linked open data modeling language"
category = "main"
optional = false
python-versions = ">=3.7.1,<4.0.0"
```

ssh -L 1111:<dbhost>:5432 -o PreferredAuthentications=password -o PubkeyAuthentication=no <user>>@<sshhost>.lbl.gov

may require some system or python resources for accessing Postgres databases

utils/pgsql2sqlite.sh: 5: [[: not found

sqlite3 resources/felix_dump.db "update auth_user set password = NULL;"
Error: in prepare, no such table: auth_user (1)

set blast_thread_count in project.Makefile

requires blastx on system path

update jinja templates
https://github.com/linkml/linkml/tree/main/linkml/generators/docgen?

resources/swiss_entries.json
writing to yaml MUCH slower than writing to json 

have to keep enum pvs up to date with data (or vice versa!)
