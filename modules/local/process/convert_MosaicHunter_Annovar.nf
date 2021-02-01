/*
* Convert MosaicHunter to Annovar input
*/
process CONVERT_MOSAICHUNTER_ANNOVAR {

  tag "$sampID"

  input:
  tuple val(sampID), file(tsv)

  output:
  tuple val(sampID), file("${sampID}.avinput.txt")

  script:
  """
awk  '{OFS = "\t"}
{if (\$3=\$7)
      print \$1, \$2, \$2, \$3, \$9, \$10/\$4, \$4, \$8"/"\$10;
else
    print \$1, \$2,  \$2, \$3, \$7,  \$8/\$4, \$4, \$10"/"\$8;
   }'    $tsv > ${sampID}.avinput.txt
  """

}
