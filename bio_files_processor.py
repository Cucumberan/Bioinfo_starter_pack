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
