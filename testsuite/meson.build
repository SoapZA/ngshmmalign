#########
# tests #
#########

foreach i : [
  'dna_array_compile_test',
  'type_caster_compile_test',
  'reference_compile_test',
  'hmmalign_forward_compile_test']
    test(
      i,
      executable(
        i + '_exe', files([
          i + '.cpp']),
        install : true,
        dependencies : ngshmmalign_deps,
        include_directories : ngshmmalign_include_dir)
      )
endforeach
