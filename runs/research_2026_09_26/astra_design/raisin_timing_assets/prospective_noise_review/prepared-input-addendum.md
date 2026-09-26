# Prepared-input and support addendum

The first draft `inputs/cadence.simlib` was already caught by root as malformed: splitting on `S:` collided with `NOBS:`. The original bytes and failed parse are retained; this review independently encountered the same failed draft and does not claim a new historical or scientific failure.

The separately written `inputs-v2` passes this review's 11-hash metadata check, 117-row/NOBS identity, six JH rows, and required header keys. Its phase span is −14.9991611 to +44.9086155 days at the declared R4-representable peak/redshift. The combined table bounds imply a same-row peak-displacement interval of [−36.4577819,+7.2662189] observer days. These are tabulation support limits, not empirical validation of the template.

The output-only all-USRFUN support patch passes static review; ordinary/instrumented exact replay and observed all-call support counts are still necessary engineering gates. No native jobs were run by this reviewer. The prospective final protocol is owned by AstraSED/root; this note neither releases execution nor certifies a final science result.
