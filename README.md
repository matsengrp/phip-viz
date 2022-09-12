# phip-viz
Visualization of PhIPseq datasets

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Docker Repository on Quay](https://quay.io/repository/matsengrp/phippery/status "Docker Repository on Quay")](https://quay.io/repository/matsengrp/phippery)
[![build and test](https://github.com/matsengrp/phippery/workflows/build%20and%20test/badge.svg)](https://github.com/matsengrp/phippery/blob/master/.github/workflows/build-and-test.yaml)

Please see the
[documentation](https://matsengrp.github.io/phippery/)
for further details.

## Docker quickstart (recommended)

```
cd pickle_data/ # a place where .phip binaries exist
docker run -p 8501:8501 -v $PWD:/app/data/ quay.io/matsengrp/phip-viz
```

