#
# rnarry.nrclip.PipelineControl
#  - Provides utility functions for process control
#
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

import os
import inspect
import string
import tempfile
import shutil

from rnarry.nrclip import Paths

__all__ = ['runproc', 'ExternalProcessError', 'TemporaryDirectory',
           'for_each_sample']

class ExternalProcessError(Exception):
    pass

def runproc(origcmd, ignore_error=False):
    callerframe = inspect.getouterframes(inspect.currentframe())[1]

    cmdtemplate = string.Template(origcmd)
    tmplvalues = Paths.__dict__.copy()
    tmplvalues.update(callerframe[0].f_globals)
    tmplvalues.update(callerframe[0].f_locals)
    command = cmdtemplate.substitute(tmplvalues).strip()

    ret = os.system(command)

    if ret != 0:
        raise ExternalProcessError(
                "[%s] Process returned %d while running command %s" %
                (callerframe[3], ret, repr(command)))


class TemporaryDirectory(object):
    def __init__(self, dir=Paths.TMP_DIR):
        self.dir = dir
        self.path = None

    def __enter__(self):
        self.path = tempfile.mkdtemp(dir=self.dir)
        return self.path

    def __exit__(self, type, value, traceback):
        if self.path is not None:
            shutil.rmtree(self.path)

def for_each_sample(inputpat, outputpat, samples):
    return [[inputpat(sample), outputpat(sample), sample]
            for sample in samples]

