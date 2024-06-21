#!/bin/sh

testwigpeNoExtOut()
{
  ${rangerCmd} wigpe ${testBam} ${testBamWigoutNoExt}  >${stdoutF} 2>${stderrF}
  rtrn=$?
  th_assertTrueWithNoOutput ${rtrn} "${stdoutF}" "${stderrF}"
  assertTrue 'wig result file has bad filename' "[ -e ${goodWig} ]"
}

testwigpeWithExtOut()
{
  ${rangerCmd} wigpe ${testBam} ${testBamWigout}  >${stdoutF} 2>${stderrF}
  rtrn=$?
  th_assertTrueWithNoOutput ${rtrn} "${stdoutF}" "${stderrF}"
  assertTrue 'wig result file has an extra .wig extension' "[ -e ${goodWig} ]"
}

testSXZ()
{
  ${rangerCmd} wigpe -sxz ${testBam} ${testBamWigout}  >${stdoutF} 2>${stderrF}
  rtrn=$?
  th_assertTrueWithNoOutput ${rtrn} "${stdoutF}" "${stderrF}"
  assertTrue 'expecting *_chr10_Pos.wig.gz ' "[ -e ${oriName}_chr10_Pos.wig.gz ]"
  assertTrue 'expecting *_chr10_Neg.wig.gz ' "[ -e ${oriName}_chr10_Neg.wig.gz ]"
}

testXZ()
{
  ${rangerCmd} wigpe -xz ${testBam} ${testBamWigout}  >${stdoutF} 2>${stderrF}
  rtrn=$?
  th_assertTrueWithNoOutput ${rtrn} "${stdoutF}" "${stderrF}"
  assertTrue 'expecting *_Pos.wig.gz ' "[ -e ${oriName}_Pos.wig.gz ]"
  assertTrue 'expecting *_Neg.wig.gz ' "[ -e ${oriName}_Neg.wig.gz ]"
}

testSZ()
{
  ${rangerCmd} wigpe -sz ${testBam} ${testBamWigout}  >${stdoutF} 2>${stderrF}
  rtrn=$?
  th_assertTrueWithNoOutput ${rtrn} "${stdoutF}" "${stderrF}"
  assertTrue 'expecting *_chr10.wig.gz ' "[ -e ${oriName}_chr10.wig.gz ]"
}

testSX()
{
  ${rangerCmd} wigpe -xs ${testBam} ${testBamWigout}  >${stdoutF} 2>${stderrF}
  rtrn=$?
  th_assertTrueWithNoOutput ${rtrn} "${stdoutF}" "${stderrF}"
  assertTrue 'expecting *_chr10_Pos.wig ' "[ -e ${oriName}_chr10_Pos.wig ]"
  assertTrue 'expecting *_chr10_Neg.wig ' "[ -e ${oriName}_chr10_Neg.wig ]"
}

testZ()
{
  ${rangerCmd} wigpe -z ${testBam} ${testBamWigout}  >${stdoutF} 2>${stderrF}
  rtrn=$?
  th_assertTrueWithNoOutput ${rtrn} "${stdoutF}" "${stderrF}"
  assertTrue 'expecting *.wig.gz ' "[ -e ${oriName}.wig.gz ]"
}

testX()
{
  ${rangerCmd} wigpe -x ${testBam} ${testBamWigout}  >${stdoutF} 2>${stderrF}
  rtrn=$?
  th_assertTrueWithNoOutput ${rtrn} "${stdoutF}" "${stderrF}"
  assertTrue 'expecting *_Pos.wig ' "[ -e ${oriName}_Pos.wig ]"
  assertTrue 'expecting *_Neg.wig ' "[ -e ${oriName}_Neg.wig ]"
}

testS()
{
  ${rangerCmd} wigpe -s ${testBam} ${testBamWigout}  >${stdoutF} 2>${stderrF}
  rtrn=$?
  th_assertTrueWithNoOutput ${rtrn} "${stdoutF}" "${stderrF}"
  assertTrue 'expecting *_chr10.wig ' "[ -e ${oriName}_chr10.wig ]"
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

  rangerCmd='../bin/peakranger'  # save command name in variable to make future changes easy
  testDir="${SHUNIT_TMPDIR}/some_test_dir"
  testBam='test_data/bit.sorted.bam'
  testBamWigout='ga.wig'
  goodWig='ga.wig'
  testBamWigoutNoExt='ga'
  oriName='ga'
}

tearDown()
{
  rm -f ${oriName}*.wig
  rm -f ${oriName}*.wig*
}


# load and run shUnit2
[ -n "${ZSH_VERSION:-}" ] && SHUNIT_PARENT=$0
. shunit2/shunit2

