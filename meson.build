project(
    'nthits',
    'cpp',
    version: '0.0.1',
    license: 'MIT',
    default_options: ['cpp_std=c++17', 'warning_level=3']
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

executable(
    'nthits',
    ['src/main.cpp'],
    dependencies: dependencies,
    include_directories: include_dirs,
    install: true
)

executable(
    'nthits-bfq',
    ['src/bf_query.cpp'],
    dependencies: dependencies,
    include_directories: include_dirs,
    install: true
)