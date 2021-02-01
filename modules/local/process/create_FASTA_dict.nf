/*
 * Create fasta dictionary and index
 */

 process CREATE_FASTA_DICT {

   input:
   tuple file("ref.fasta"), file("ref.fasta.fai")

   output:
   tuple file("ref.fasta"), file("ref.fasta.fai"), file("ref.dict")

   script:
   """
   gatk CreateSequenceDictionary -R ref.fasta
   """
 }
