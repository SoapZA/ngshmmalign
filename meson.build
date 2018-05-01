project(
  'ngshmmalign',
  'cpp',
  version : '0.1.1',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++11',
    'b_ndebug=if-release'],
  license : 'GPL-2+',
  meson_version : '>= 0.46.0')

add_global_arguments([
    '-DPACKAGE_NAME="' + meson.project_name() + '"',
    '-DPACKAGE_VERSION="' + meson.project_version() + '"'],
  language : 'cpp')

################
# dependencies #
################

# boost
ngshmmalign_deps = [
  dependency(
    'boost', modules : [
      'system',
      'program_options',
      'filesystem',
      'iostreams'],
    required : true)]

# OpenMP
if get_option('openmp')
  ngshmmalign_deps += dependency('openmp', required : true)
endif

###########
# sources #
###########

subdir('src')

#########
# tests #
#########

if get_option('tests')
  subdir('testsuite')
endif
