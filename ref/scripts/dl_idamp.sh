#!/bin/bash
#
# dl_idamp.sh
#
# Download 519families_idmapping
#
# Copyright (C) 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

aws s3 cp s3://sscwgs/IDmapping/519families_idmapping .
