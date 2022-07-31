#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2022 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2022-07-31 00:59

import varcode


def get_top_effect(effs, pU_mode=False):
    """Get top effect of a EffectCollection object.

    input is varcode.EffectCollection([])
    """

    if len(effs) == 0:
        return None
    if not pU_mode:
        return effs.top_priority_effect()
    priority_types = [
        "rRNA",
        "rRNA_pseudogene",
        "Mt_rRNA",
        "tRNA",
        "Mt_tRNA",
        "snoRNA",
        "snRNA",
        "scaRNA",
        "scRNA",
        "vault_RNA",
        "miRNA",
    ]
    priority_dict = {t: [] for t in priority_types}
    priority_dict["others"] = []
    for eff in effs:
        if (
            hasattr(eff, "gene")
            and hasattr(eff.gene, "biotype")
            and eff.gene.biotype in priority_types
        ):
            priority_dict[eff.gene.biotype].append(eff)
        elif (
            hasattr(eff, "transcript")
            and hasattr(eff.transcript, "biotype")
            and eff.transcript.biotype in priority_types
        ):
            priority_dict[eff.transcript.biotype].append(eff)
        else:
            priority_dict["others"].append(eff)
    for eff_list in priority_dict.values():
        if len(eff_list) > 0:
            return varcode.EffectCollection(eff_list).top_priority_effect()
    return None
