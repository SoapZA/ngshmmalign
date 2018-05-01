ngshmmalign_include_dir = include_directories('.')

###############
# executables #
###############

ngshmmalign = executable(
  'ngshmmalign', files([
    'main.cpp']),
  install : true,
  dependencies : ngshmmalign_deps,
  include_directories : ngshmmalign_include_dir)
