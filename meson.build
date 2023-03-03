project(
    'nthits',
    'cpp',
    version: '1.0.0',
    license: 'MIT',
    default_options: ['cpp_std=c++17', 'warning_level=3', 'werror=true']
)

configure_file(
  input : 'include/nthits.hpp',
  output : 'nthits.hpp',
  configuration : {
    'PROJECT_VERSION' : meson.project_version()
  }
)

compiler = meson.get_compiler('cpp')

dependencies = [
    dependency('openmp'),
    compiler.find_library('btllib')
]

include_dirs = [
    include_directories('include'),
    include_directories('vendor/argparse/include'),
]

tests_files = [
    'tests/min_thresh.cpp'
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

# TESTS ------------------------------------------------------------------------
catch2 = dependency('Catch2', method : 'cmake', modules : ['Catch2::Catch2WithMain'], required: false)

if not catch2.found()
  error('Tests require catch2 to be installed.')
endif

foreach file : tests_files
    t = file.split('.cpp')[0].split('/')[-1]
    test(t, executable(t + '_test', file,
        include_directories : include_dirs,
        dependencies : dependencies + [catch2])
    )
endforeach
