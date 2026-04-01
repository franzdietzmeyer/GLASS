"""
Unit tests for ptm_terminal_filter.filter_ptm_terminal_by_chain_order (Rosetta PTM termini).
"""
import pytest

from ptm_terminal_filter import filter_ptm_terminal_by_chain_order


class TestFilterPtmTerminalByChainOrder:
    """Per-chain ordinal filter: first/last 4 residues excluded (default)."""

    def test_short_chain_returns_empty(self):
        # Need L >= 9 for exclude_each_end=4
        pdb_per = list(range(1, 9))  # L=8
        assert filter_ptm_terminal_by_chain_order(pdb_per, [5], exclude_each_end=4) == []

    def test_min_length_chain_boundary(self):
        pdb_per = list(range(1, 10))  # L=9, valid ordinals 5 only (first_ord=5, last_ord=5)
        assert filter_ptm_terminal_by_chain_order(pdb_per, [5], exclude_each_end=4) == [5]
        assert filter_ptm_terminal_by_chain_order(pdb_per, [1, 9], exclude_each_end=4) == []

    def test_twenty_residue_chain_excludes_terminal_candidates(self):
        # PDB nums 1..20, ordinal = index; allowed ordinals 5..16
        pdb_per = list(range(1, 21))
        got = filter_ptm_terminal_by_chain_order(pdb_per, [1, 4, 5, 10, 16, 17, 20], exclude_each_end=4)
        assert got == [5, 10, 16]

    def test_two_chains_independent_overlapping_pdb_numbers(self):
        """Same PDB labels on two chains: filtering uses each chain's list only (not global min/max)."""
        chain_a = list(range(10, 30))  # 20 residues, PDB 10..29
        chain_b = list(range(10, 30))  # identical labels, different molecule
        a_keep = filter_ptm_terminal_by_chain_order(chain_a, [10, 14, 25, 29], exclude_each_end=4)
        b_keep = filter_ptm_terminal_by_chain_order(chain_b, [10, 14, 25, 29], exclude_each_end=4)
        # Ordinals 1-4: PDB 10-13 (drop 10); ordinals 17-20: PDB 26-29 (drop 29). PDB 25 is ordinal 16 (keep).
        assert a_keep == [14, 25]
        assert b_keep == [14, 25]

    def test_unknown_candidate_dropped(self):
        pdb_per = list(range(1, 21))
        assert filter_ptm_terminal_by_chain_order(pdb_per, [999], exclude_each_end=4) == []

    def test_custom_exclude_each_end(self):
        pdb_per = list(range(1, 11))  # L=10, exclude 2 each end -> first_ord=3, last_ord=8
        got = filter_ptm_terminal_by_chain_order(pdb_per, [1, 2, 5, 9, 10], exclude_each_end=2)
        assert got == [5]
