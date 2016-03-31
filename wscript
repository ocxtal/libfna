#! /usr/bin/env python
# encoding: utf-8

def options(opt):
	opt.load('compiler_c')

def configure(conf):
	conf.recurse('zf')

	conf.load('ar')
	conf.load('compiler_c')
	conf.env.append_value('CFLAGS', '-O3')
	conf.env.append_value('CFLAGS', '-march=native')


def build(bld):
	bld.recurse('zf')

	bld.stlib(
		source = ['fna.c'],
		target = 'fna',
		lib = bld.env.LIB_ZF,
		use = ['zf'])

	bld.program(
		source = ['fna.c'],
		target = 'unittest',
		lib = bld.env.LIB_ZF,
		use = ['zf'],
		defines = ['TEST'])
