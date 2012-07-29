#!/usr/bin/env python
#
# Copyright (C) 2012 Hyeshik Chang
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#

from ruffus import *
from rnarry.nrclip import (
    Paths, DataPreparation, SequenceProcessing, ContaminantFilter,
    SequenceAlignment, SequenceAnnotation, SequenceMasking,
    DerivedDatabaseBuilding, TranscriptomeAnalysis, GeneralStatistics)
from rnarry.nrclip import Options
from itertools import chain
import os

task_modules = [
    DataPreparation,
    SequenceProcessing,
    ContaminantFilter,
    SequenceAlignment,
    SequenceAnnotation,
    SequenceMasking,
    DerivedDatabaseBuilding,
    TranscriptomeAnalysis,
    GeneralStatistics
]

all_tasks = list(chain(*[mod.tasks() for mod in task_modules]))

# Simplify module names in diagrams
for task in all_tasks:
    module = task.__module__.split('.')[-1]
    task.pipeline_task.display_name = '%s\\n%s' % (module, task.__name__)

# Prepare working directories
for dirpath in Paths.ALL_SUBDIRS:
    if not os.path.isdir(dirpath):
        os.mkdir(dirpath)

pipeline_printout_graph('flowchart.pdf', 'pdf', all_tasks,
                        pipeline_name='Pipeline for LIN28A CLIP-seq and RPF')
#pipeline_printout_graph('flowchart.dot', 'dot', all_tasks,
#                        pipeline_name='Pipeline for LIN28A CLIP-seq and RPF')
pipeline_run(all_tasks, verbose=5, multiprocess=Options.NUM_PARALLEL)

# ex: ts=8 sts=4 sw=4 et
