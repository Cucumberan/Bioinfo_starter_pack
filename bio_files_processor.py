def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None) -> None:
    """
    Converts a multiline FASTA file to a one-line FASTA file.

    Args:
        input_fasta (str): The path to the input multiline FASTA file.
        output_fasta (str, optional): The path to the output one-line FASTA file.
            If not provided, the output file will have the same name as the input file
            with '_one_line.fasta' appended.

    Returns:
        None
    """
    if not output_fasta:
        output_fasta = input_fasta.split('.')[0] + '_one_line.fasta'

    with open(input_fasta, 'r') as multiline_fasta, open(output_fasta, 'w') as one_line_fasta:
        current_sequence = ''

        for line in multiline_fasta:
            if line.startswith('>'):
                if current_sequence:
                    one_line_fasta.write(current_sequence + '\n')
                one_line_fasta.write(line)
                current_sequence = ''
            else:
                current_sequence += line.strip()
        if current_sequence:
            one_line_fasta.write(current_sequence + '\n')


def select_genes_from_gbk_to_fasta(input_gbk: str, genes: list, output_fasta: None) -> None:
    """
    Extracts gene information from a GenBank file and saves it in FASTA format. For now, this function finds only
    those protein sequences that correspond to genes that are taken as input in the form of a list of genes (gene
    symbols). But someday it may become a full-fledged gdk parser. Possibly. Not sure, but still hope.

        Args:
            input_gbk (str): Path to the input GenBank file.
            genes (list): List of gene names to extract.
            output_fasta (str, optional): Path to the output FASTA file. If not provided,
                a default name will be generated based on the input GenBank file name.
        Returns:
            None

        """
    start_substring = 'translation='
    end_substring = '"'
    if not output_fasta:
        output_fasta = input_gbk.split('.')[0] + '.fasta'

    with open(input_gbk, 'r') as gbk, open(output_fasta, 'w') as output_file:
        found_gene = False
        found_start = False
        found_end = False
        for gene in genes:
            output_file.write(gene + '\n')
            for line in gbk:
                if gene in line:
                    found_gene = True
                if found_gene:
                    if found_start or start_substring in line:
                        found_start = True
                        if found_end:
                            if end_substring in line:
                                output_file.write(line.strip() + '\n')
                                break
                            output_file.write(line.strip() + '\n')
                        elif end_substring in line:
                            found_end = True
                            output_file.write(line.strip() + '\n')
            found_gene = False
            found_start = False
            found_end = False
            gbk.seek(0)
            output_file.write('\n')
