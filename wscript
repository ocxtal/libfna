#! /usr/bin/env python
# encoding: utf-8

def options(opt):
	opt.load('compiler_c')

def configure(conf):
	conf.load('ar')
	conf.load('compiler_c')
	if 'LIB_Z' not in conf.env:
		conf.check_cc(lib = 'z')

	conf.env.append_value('CFLAGS', '-O3')
	conf.env.append_value('CFLAGS', '-march=native')


def build(bld):
	bld.stlib(
		source = ['fna.c', 'kopen.c'],
		target = 'fna',
		lib = ['z'])

	bld.program(
		source = ['fna.c', 'kopen.c'],
		target = 'unittest',
		lib = ['z'],
		defines = ['TEST'])
