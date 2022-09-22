project(
    'nthits',
    'cpp',
    version: '0.0.1',
    license: 'MIT',
    default_options: ['cpp_std=c++17']
)

compiler = meson.get_compiler('cpp')

dependencies = [
    dependency('openmp'),
    compiler.find_library('btllib')
]

include_dirs = [
    include_directories('include'),
    include_directories('vendor/argparse/include')
]

sources = run_command('scripts/list_sources.sh', check: true).stdout().strip().split()

executable(
    'nthits',
    ['src/main.cpp'] + sources,
    dependencies: dependencies,
    include_directories: include_dirs,
    install: true
)