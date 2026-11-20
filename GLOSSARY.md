# Rasusa

Rasusa randomly subsamples sequencing data. It answers one question — *which subset of
this data should I keep?* — for two kinds of input: unaligned reads and aligned
sequences.

## Language

### Data

**Read**:
A single sequence produced by a sequencer, with its identifier and (usually) per-base
quality scores.
_Avoid_: Sequence, record, entry

**Template**:
The unit of selection. A template is one read for single-end data, or a read and its
mate for paired data. Both mates of a template are always kept or dropped together —
never split.
_Avoid_: Pair, fragment, cluster

**Name-grouped**:
A property of a file in which all reads belonging to the same template appear
consecutively. Sorting *by* name implies it, but a file can be name-grouped without
being sorted.
_Avoid_: Name-sorted, collated, query-sorted

**Mapped**:
A read that has been assigned a position on a reference sequence. Subsampling mapped
data is a distinct problem from subsampling reads, and Rasusa treats them separately.
_Avoid_: Aligned, placed

### Selection

**Subsample**:
To select a random subset of templates from the input and discard the rest. Rasusa never
modifies, trims, or reorders what it keeps.
_Avoid_: Downsample, filter, sample, thin

**Target**:
What the user asked for — a depth of coverage, a total number of bases, a number of
reads, or a fraction of the input. Targets say *how much* data to keep, and are
independent of how Rasusa goes about it.
_Avoid_: Goal, amount, quantity

**Strategy**:
How Rasusa reaches a target — how many times it reads the input and what it must know
before it can start writing. Orthogonal to the target: the same target may be reachable
by more than one strategy, and the choice between them is a speed/accuracy trade.
_Avoid_: Mode, method, approach

**Two-pass**:
The strategy that reads the input twice: once to measure it, then again to write the
selection. Measuring first is what makes exactness possible.
_Avoid_: Double-pass, exact mode

**One-pass**:
The strategy that reads the input once, deciding each template's fate as it arrives.
Cannot know the size of the input, and so cannot hit a target exactly.
_Avoid_: Streaming mode, single-pass, online

**Exact**:
A guarantee that the output matches the target precisely — the requested number of
reads, or the requested fraction of them, with no error. Only two-pass can be exact.
_Avoid_: Precise, deterministic (a result can be reproducible without being exact)

**Approximate**:
A result that is correct in expectation but not guaranteed. An approximate fraction
converges on the target as the input grows, so the error is largest on small inputs and
negligible on large ones.
_Avoid_: Inexact, fuzzy, lossy

### Measures

**Coverage**:
How many times over the target genome is represented by a set of reads — total bases
divided by genome size. Expressed as a multiple, e.g. `20x`.
_Avoid_: Depth, depth of coverage, redundancy

**Genome size**:
The size of the genome that coverage is calculated with respect to. Supplied by the
user; Rasusa cannot infer it from the reads.
_Avoid_: Reference size, target size

**Fraction**:
The proportion of the input to keep, as a value between 0 and 1.
_Avoid_: Percentage, proportion, ratio, rate

**Seed**:
The value that fixes Rasusa's random choices, so that the same input, target, strategy
and seed always produce the same output. Reproducibility is a property of the seed and
is independent of whether the result is exact.
_Avoid_: Random state

### Behaviour

**Strict**:
The stance that an unreachable target is an error rather than a warning. Strictness is
only meaningful when a guarantee exists to enforce, so it applies to exact results only.
_Avoid_: Fail-fast, pedantic mode
