# Copyright (c) 2015 Boocock James <james.boocock@otago.ac.nz>
# Author: Boocock James <james.boocock@otago.ac.nz>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import os
import subprocess
import shlex
import logging

from pyrallel import *
from fine_mapping_pipeline.expections.error_codes import *
logger = logging.getLogger(__name__)

def run_command_return_output(command,cleanup=None ,error=GENERIC_PROCESS_FAILURE, shell=False):
    """ 
        Runs one command and returns the stdout to the user.
    """
    devnull = open(os.devnull, 'w')
    logging.info("Shell {0}".format(shell))
    try:
        logging.info('Running command: {}'.format(command))
        if not shell:
            command = shlex.split(command)
            output = subprocess.check_output(command, stderr=devnull, shell=shell)
        else:
            output = subprocess.check_output(command, stderr=devnull, shell=shell)

    except subprocess.CalledProcessError:
        logger.error("Command: {0} failed".format(' '.join(command)))
        sys.exit(error)
    return output.decode("utf8")


def run_command(command,cleanup=None,error=GENERIC_PROCESS_FAILURE,exit_on_failure=True, stdout=None, shell=False):
    """
        Runs a command and returns the stderr to the user.
    """
    devnull = open(os.devnull, 'w')
    if stdout is None:
        stdout = devnull
    try:
        logging.info('Running command: {}'.format(command))
        if shell:
        #    stdout = open("/Users/smilefreak/test.txt", "w")
            subprocess.check_call(command,stderr=devnull, stdout=stdout, shell=True)
        else:
            command = shlex.split(command)
        #    stdout = open("/Users/smilefreak/test.txt", "w")
            subprocess.check_call(command,stderr=devnull, stdout=stdout)
    except subprocess.CalledProcessError:
        logger.error("Command: {0} failed".format(' '.join(command)))
        if exit_on_failure:
            sys.exit(error)
        logger.warning("Did not exit because of a user request")
        return False
    return True

def run_commands(commands, tool_name="" ,cores=6, stdouts=None):
    """
        Runs one command and returns the stdout the user.

        This utilises the pyrallel package.
    """   
    queue_jobs(commands, tool_name, threads=6,stdouts=None)
