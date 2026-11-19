---
status: accepted
---

# One-pass subsampling is fraction-only

Rasusa's two-pass strategy measures the input before writing, which is what lets it hit a
target exactly but costs a full extra read of the file — painful on very large inputs. We
are adding a one-pass strategy, but **only for fraction targets**, implemented as
independent Bernoulli trials per template. Number-of-reads, bases, and coverage targets
remain two-pass permanently. In particular we considered and rejected reservoir sampling,
which would have made one-pass number-of-reads *exact*, because it costs memory
proportional to the number of reads kept.

## Considered options

**Bernoulli sampling for fraction targets** (chosen). Keep each template with probability
*p*. Constant memory, output streams in input order, no knowledge of the input size
required. The result is approximate — the count is binomially distributed — but the
relative error is `sqrt((1-p)/(n·p))`, which shrinks as the input grows. On a ~700 GB
FASTQ that is on the order of 0.01%; it is only material on small inputs, where two-pass
is already fast. **The error vanishes precisely where the speedup matters**, which is what
makes the trade acceptable.

**Reservoir sampling for number-of-reads targets** (rejected). Algorithm L yields exactly
*k* templates, uniformly, in a single pass — genuinely exact and genuinely one-pass, so it
is the obvious thing for a future reader to suggest. It was rejected on memory. A
reservoir must hold *k* whole records, not *k* indices: 500,000 Illumina reads is roughly
200 MB, and 500,000 long reads roughly 10 GB, against the ~1.25 MB the two-pass strategy
uses for ten million reads. It also cannot emit anything until end-of-input, since any
later template may evict an earlier one, making it a poor pipeline citizen. Keeping
Rasusa's memory footprint small is a stronger commitment than making this particular
target one-pass.

**Bases and coverage targets** (not possible). Both require the input's total base count
to compute a rate, which is exactly what a single pass cannot know. The streaming
alternative — keep templates until the base budget is spent — is biased toward the start
of the file, which is unacceptable for data whose quality drifts over a run.

## Consequences

- One-pass and strictness are mutually exclusive. Strictness enforces a guarantee, and
  one-pass makes none, so the combination is rejected rather than silently ignored.
- Paired input is processed in lockstep, both sides advancing together on a single random
  draw per template, because there is no precomputed selection to share between two
  separate passes.
- For unaligned alignment input, templates can only be grouped by comparing against the
  immediately preceding record, so one-pass requires name-grouped input and must guard
  against being handed something else.
- Standard input becomes achievable for fraction targets, since nothing needs to be
  re-read. That is deliberately left to later work.
- Choosing memory over exactness here is a precedent: future streaming features should
  assume constant memory is the constraint that binds.
