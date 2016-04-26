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

	conf.env.append_value('LIB_FNA', conf.env.LIB_ZF)
	conf.env.append_value('DEFINES_FNA', conf.env.DEFINES_ZF)
	conf.env.append_value('OBJ_FNA', ['fna.o'] + conf.env.OBJ_ZF)


def build(bld):
	bld.recurse('zf')

	bld.objects(source = 'fna.c', target = 'fna.o')

	bld.stlib(
		source = ['unittest.c'],
		target = 'fna',
		use = bld.env.OBJ_FNA,
		lib = bld.env.LIB_FNA,
		defines = bld.env.DEFINES_FNA)

	bld.program(
		source = ['unittest.c'],
		target = 'unittest',
		use = bld.env.OBJ_FNA,
		lib = bld.env.LIB_FNA,
		defines = ['TEST'] + bld.env.DEFINES_FNA)
