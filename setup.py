#!/usr/bin/env python
"""Description
Setup script for "ncHMR detector: a computational framework to systematically reveal non-classical functions of histone-modification regulators"
Copyright (c) 2019 Shengen Hu <tarelahu@gmail.com>
This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).
"""
import os
import sys
import subprocess
import platform
from distutils.core import setup#, Extension
import distutils.command.install_lib
#from setuptools import setup, find_packages


def sp(cmd):
    '''
    Call shell cmd or software and return its stdout
    '''
    a=subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell='TRUE')
    ac = a.communicate()
    return ac

def compile_bedtools():
    curdir = os.getcwd()
    os.chdir('refpackage/bedtools')
    sp('make 1>/dev/null 2>&1 ')
    #sp('chmod 755 *')
    os.chdir(curdir)
    
def check_bedtools():
    checkhandle = sp('which bedtools')
    if checkhandle[0].strip() == "":
        return 0
    else:
        return 1
def check_R():
    checkhandle = sp('which Rscript')
    if checkhandle[0].strip() == "":
        return 0
    else:
        return 1   


class my_install_lib(distutils.command.install_lib.install_lib):
    def run(self):
        distutils.command.install_lib.install_lib.run(self)
        mode = 755
        # here we start with doing our overriding and private magic ..
        for filepath in self.get_outputs():
            if "bigWigSummary" in filepath or "bedtools" in filepath:
            #if self.install_scripts in filepath:
            #    log.info("Overriding setuptools mode of scripts ...")
            #    log.info("Changing ownership of %s to uid:%s gid %s" %
            #             (filepath, uid, gid))
            #    os.chown(filepath, uid, gid)
            #    log.info("Changing permissions of %s to %s" %
            #             (filepath, oct(mode)))
                os.chmod(filepath, mode)

def main(): 
#    if sys.version_info[0] != 2 or sys.version_info[1] < 7:
#	    print >> sys.stderr, "ERROR: ncHMR_detector requires Python 2.7"
#	    sys.exit()
    has_R = check_R()
    if has_R == 0:
	    print("ERROR: ncHMR_detector requires R & Rscript under default PATH", file=sys.stderr)
	    sys.exit()

    OS = platform.system()
    if OS == "Linux":
        bwsum_software = "bigWigSummary_linux"
    elif OS == "Darwin":
        bwsum_software = "bigWigSummary_mac"
    else:
        wlog("detected system is nither linux nor mac, try linux version of bigWigSummary",logfile)
        bwsum_software = "bigWigSummary_linux"
        
#    has_bedtools = check_bedtools()
    has_bedtools=0
    print('Intalling ncHMR_detector, may take "several" minutes')
    if has_bedtools == 0:
        compile_bedtools()
        sp('mv refpackage/bedtools/bin/bedtools lib/')
        setup(name="HMRpipe",
              version="1.3",
              description="ncHMR detector: a computational framework to systematically reveal non-classical functions of histone-modification regulators",
              author='Shengen Hu',
              author_email='Tarelahu@gmail.com',
              url='https://github.com/Tarela/ncHMR_detector.git',
              package_dir={'HMRpipe' : 'lib'},
              packages=['HMRpipe'],
              package_data={'HMRpipe': ['%s'%bwsum_software,'bedtools',#'Config/template.conf',
                                      #'Rscript/analysis.r',
                                      #'Rscript/individual_qc.r',
                                      #'Rscript/readsbulkQC.r',
                                      #'Rscript/detectNonCanonical.r'
                                         ]},
              scripts=['bin/ncHMR_detector_py3'],#,'refpackage/bedtools/bin/bedtools','refpackage/bwsummary/%s'%bwsum_software],
                        
              classifiers=[
            'Development Status :: version1.3 finish',
            'Environment :: Console',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Artistic License',
            'Operating System :: POSIX',
            'Programming Language :: Python',
            'Topic :: pipeline',
            ],
              requires=[],
              cmdclass={'install_lib':my_install_lib}
          )
        print('bedtools is not detected under default PATH, bedtools is also installed')
        #print('Installation of ncHMR_detector is DONE')
    
    else:
        setup(name="HMRpipe",
              version="1.3",
              description="ncHMR detector: a computational framework to systematically reveal non-classical functions of histone-modification regulators",
              author='Shengen Hu',
              author_email='Tarelahu@gmail.com',
              url='https://github.com/Tarela/ncHMR_detector.git',
              package_dir={'HMRpipe' : 'lib'},
              packages=['HMRpipe'],
              #package_data={}
              package_data={'HMRpipe': ['%s'%bwsum_software]},#,#'Config/template.conf',
                                      #'Rscript/analysis.r',
                                      #'Rscript/individual_qc.r',
                                      #'Rscript/readsbulkQC.r',
                                      #'Rscript/detectNonCanonical.r'
                                      #   ]},
              #scripts=['bin/ncHMR_detector_py3','refpackage/bwsummary/%s'%bwsum_software],
              scripts=['bin/ncHMR_detector_py3'],#,'refpackage/bwsummary/%s'%bwsum_software],
              #data_files=[('/Users/sh8tv/bin',['refpackage/bwsummary/%s'%bwsum_software])],
              classifiers=[
            'Development Status :: version1.3 finish',
            'Environment :: Console',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Artistic License',
            'Operating System :: POSIX',
            'Programming Language :: Python',
            'Topic :: pipeline',
            ],
              requires=[],
              cmdclass={'install_lib':my_install_lib}
          )

    #import HMRpipe

    print('Installation of ncHMR_detector is DONE')


if __name__ == '__main__':
    main()



