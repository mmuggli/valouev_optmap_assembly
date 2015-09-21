# Notes by Martin D. Muggli

This repository is for maintanance/enhancements to the software accompanying this publication:

Anton Valouev, David C. Schwartz, Shiguo Zhou, and Michael S. Waterman. "An algorithm for assembly of ordered restriction maps from single DNA molecules"  October 24, 2006.  PNAS. no. 43 www.pnas.org/cgi/doi/10.1073/pnas.0604040103 vol. 103   

# Usage
This tool doesn't quite take the output of the aligner.  Two tools in the contrib directory help fill in the gaps.

Compile with g++ assmbly.cc

    $ valuev_optmap_assembly/contrib/gen_skipped.py sim_ecoli_XhoI_rev.valuev sim_ecoli_XhoI_rev.ovlps sim_ecoli_XhoI_rev.skipped
    $ valuev_optmap_assembly/contrib/patchin_hash_keys.py sim_ecoli_XhoI_rev.valuev sim_ecoli_XhoI_rev.ovlps sim_ecoli_XhoI_rev.ovlps.patched
    $ valuev_optmap_assembly/a.out 

n.b. You must modify strings in assmbly.cc to match your input/output files. For this example:

    map_filename =  "sim_ecoli_XhoI_rev.valuev"; // the reads
    ovlp_filename =   "sim_ecoli_XhoI_rev.ovlps.patched"; // the alignments with read position in the read file prefixed
    comp_prefix = "sim_ecoli_XhoI_rev.contigs"; // the output filename prefix,
    skipped_maps_file_name =  "sim_ecoli_XhoI_rev.skipped"; // reads not found in any alignment, assumed to be chimera

I was able to simulate a set of ecoli read, align it with the overlap [aligner](https://github.com/mmuggli/valuev_optmap_alignment) "ovlp", assemble it with this tool, and then align the assembly back to an *in silico* digest of the reference genome with the "fit" aligner in the aligner repo.
