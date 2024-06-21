#!/bin/sh

test_stable_peak()
{
  for i in {1..100}
  do
    #chr11 87067606  87067845  5.18579e-05 +
    opk=87067606
    ${rangerCmd} ${format} ${data} ${input} ${out} 1>${stdoutF} 2>${stderrF}
    pk=`grep ${opk} ${out}_region.bed | cut -f 2`
    assertEquals 'Peak ${opk} should always exist' "${opk}" "${pk}"
  done
}

#-----------------------------------------------------------------------------
# suite functions
#

th_assertTrueWithNoOutput()
{
  th_return_=$1
  th_stdout_=$2
  th_stderr_=$3

  assertFalse 'unexpected output to STDOUT' "[ -s '${th_stdout_}' ]"
  assertFalse 'unexpected output to STDERR' "[ -s '${th_stderr_}' ]"

  unset th_return_ th_stdout_ th_stderr_
}

oneTimeSetUp()
{
  outputDir="${SHUNIT_TMPDIR}/output"
  mkdir "${outputDir}"
  stdoutF="${outputDir}/stdout"
  stderrF="${outputDir}/stderr"

  rangerCmd='../bin/peakranger ranger'  
  testDir="${SHUNIT_TMPDIR}/some_test_dir"
  data='test_data/chr11.region.bam'
  format='--format bam'
  input='test_data/chr11.region.input.bam'
  out='ga'
}

tearDown()
{
  rm -f ${out}_*
}


# load and run shUnit2
[ -n "${ZSH_VERSION:-}" ] && SHUNIT_PARENT=$0
. shunit2/shunit2

