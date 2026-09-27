"""Typed numerical-target descriptor; no observational qualifier implemented.

This separates an accuracy2 native configuration from the original sampled
proposal identity. Constructing the descriptor alone never admits a posterior.
The design specifies the later full2000-record qualifier; old accuracy1
consumers must not infer this target from a renamed file.
"""
from dataclasses import dataclass
import json

from native_accuracy_correction import assert_only_boost_change,identity
from target_identity import canonical


@dataclass(frozen=True)
class NativeAccuracyTarget:
    numerical_accuracy: int
    target_identity: str
    parent_proposal_identity: str
    configuration_json: str
    lineage_json: str
    scope: str = 'configuration_descriptor_only_not_posterior_qualification'

    @property
    def configuration(self):
        return json.loads(self.configuration_json)

    @property
    def lineage(self):
        return json.loads(self.lineage_json)


def target_descriptor(native1,native2,parent_proposal_identity,asset_sha256,versions,source_sha256):
    assert_only_boost_change(native1,native2)
    assert parent_proposal_identity and asset_sha256 and versions and source_sha256
    config=canonical(native2)
    lineage={'schema':'native-accuracy2-target-descriptor-v1',
             'parent_proposal_identity':parent_proposal_identity,
             'native_accuracy1_configuration':canonical(native1),
             'native_accuracy2_configuration':config,
             'asset_sha256':dict(asset_sha256),'versions':dict(versions),
             'source_sha256':dict(source_sha256),'numerical_accuracy':2,
             'scope':'configuration_descriptor_only_not_posterior_qualification'}
    return NativeAccuracyTarget(2,identity(lineage),parent_proposal_identity,
                                json.dumps(config,sort_keys=True,allow_nan=False),
                                json.dumps(lineage,sort_keys=True,allow_nan=False))
