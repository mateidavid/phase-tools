SHELL := /bin/bash

INSTALL_DIR = ../../build

.PHONY: all list clean test.add-parent-phasing test.extend-phase-set test.ngs-phase
.SECONDARY:
.DELETE_ON_ERROR:

TARGETS = test.add-parent-phasing test.extend-phase-set test.ngs-phase

all: ${TARGETS}
	@echo "all good!"

list:
	@echo "TARGETS=${TARGETS}"

clean:
	rm -f test.ngs-phase.0.in.bam test.ngs-phase.0.in.bam.bai test.ngs-phase.0.in.fa.fai test.ngs-phase.0.real-out.vcf

clean-all: clean

test.add-parent-phasing: ${INSTALL_DIR}/bin/add-parent-phasing
	diff <(${INSTALL_DIR}/bin/add-parent-phasing --child HG002 --parent HG003 HG004 -i test.add-parent-phasing.0.in.vcf 2>/dev/null) test.add-parent-phasing.0.out.vcf
	@echo "$@: all tests passed"

test.extend-phase-set: ${INSTALL_DIR}/bin/extend-phase-set
	diff <(${INSTALL_DIR}/bin/extend-phase-set -i test.extend-phase-set.0.in.vcf --base-id HG002_base --extra-id HG002_extra 2>/dev/null) test.extend-phase-set.0.out.vcf
	diff <(${INSTALL_DIR}/bin/extend-phase-set -i test.extend-phase-set.1.in.vcf --base-id HG002 --extra-id HG003 --extra-is-parent 0 2>/dev/null) test.extend-phase-set.1.out.vcf
	diff <(${INSTALL_DIR}/bin/extend-phase-set -i test.extend-phase-set.2.in.vcf --base-id HG002 --extra-id HG003 --extra-is-parent 0 2>/dev/null) test.extend-phase-set.2.out.vcf
	@echo "$@: all tests passed"

%.bam.bai: %.sam
	samtools view $< -b >$(@:.bai=)
	samtools index $(@:.bai=)

%.fa.fai: %.fa
	samtools faidx $<

test.ngs-phase: ${INSTALL_DIR}/bin/ngs-phase test.ngs-phase.0.in.bam.bai test.ngs-phase.0.in.fa.fai
	${INSTALL_DIR}/bin/ngs-phase --ref test.ngs-phase.0.in.fa --var test.ngs-phase.0.in.vcf --map test.ngs-phase.0.in.bam --sample HG002 --out test.ngs-phase.0.real-out.vcf
	diff test.ngs-phase.0.real-out.vcf test.ngs-phase.0.out.vcf
	@echo "$@: all tests passed"
