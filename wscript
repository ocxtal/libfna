#! /usr/bin/env python
# encoding: utf-8

def options(opt):
	opt.recurse('zf')
	opt.load('compiler_c')

def configure(conf):
	conf.recurse('zf')

	conf.load('ar')
	conf.load('compiler_c')
	conf.env.append_value('CFLAGS', '-O3')
	conf.env.append_value('CFLAGS', '-std=c99')
	conf.env.append_value('CFLAGS', '-march=native')


def build(bld):
	bld.recurse('zf')

	bld.stlib(
		source = ['fna.c'],
		target = 'fna',
		lib = bld.env.LIB_ZF,
		use = ['zf'])

	bld.program(
		source = ['unittest.c'],
		target = 'unittest',
		linkflags = ['-all_load'],
		use = ['fna', 'zf'],
		lib = bld.env.LIB_ZF,
		defines = ['TEST'] + bld.env.DEFINES_ZF)
