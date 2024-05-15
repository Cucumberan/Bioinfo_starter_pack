#HW_6
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

#HW_15
class OpenFasta:
    """
    A context manager class for opening and iterating through records in a FASTA file.

    This class provides a convenient way to handle the operations on a FASTA file, each
    record in a FASTA file consists of a header line that starts with '>', followed by one
    or more lines of sequence data.

    Attributes:
        file_path (str): Path to the FASTA file.
        mode (str): Mode in which the file is opened, default is 'r' for read-only.
        file (file): The file object, initialized as None and set when entering the context.

    Methods:
        __enter__(): Opens the file and returns the instance itself as a context manager.
        __exit__(exc_type, exc_value, traceback): Closes the file upon exiting the context.
        __iter__(): Returns self, making this class an iterator.
        __next__(): Returns the next sequence record from the file as a FastaRecord object.
        read_record(): Returns the next record in the file, or None if no more records.
        read_records(): Returns a list of all records in the file.

    Usage example:
        with OpenFasta('path/to/file.fasta') as fasta_file:
            for record in fasta_file:
                print(record)
    """

    def __init__(self, file_path, mode='r'):
        self.file_path = file_path
        self.mode = mode
        self.file = None

    def __enter__(self):
        self.file = open(self.file_path, self.mode)      
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.file:
            self.file.close()
    
    def __iter__(self):
        return self
    
    def __next__(self):
        if not self.file:
            raise StopIteration
        
        line = ''
        while True:
            line = self.file.readline()
            if not line:
                raise StopIteration  # The end of file
            if line.startswith('>'):
                break
        header = line.strip().split(maxsplit=1)
        id = header[0][1:]
        description = header[1] if len(header) > 1 else ''

        seq_lines = []
        while True:
            line = self.file.readline().strip()
            if not line or line.startswith('>'):
                self.file.seek(self.file.tell() - len(line) - 1)  # Return a pointer for the next __next__ call
                break
            seq_lines.append(line)
        seq = ''.join(seq_lines)
        return FastaRecord(id, description, seq)

    def read_record(self):
        try:
            return next(self)
        except StopIteration:
            return None

    def read_records(self):
        return list(self)


class FastaRecord():
    """
    Represents a single record from a FASTA file.

    This class is designed to encapsulate the information for a single FASTA record,
    which typically includes an identifier, an optional description, and the sequence itself.
    The identifier and description are parsed from the header line of the record, and the
    sequence is composed of all subsequent lines until the next header line or end of file.

    Attributes:
        id (str): The identifier of the sequence, usually a unique value.
        description (str): An optional description providing more details about the sequence.
        seq (str): The actual biological sequence data (e.g., DNA, RNA, protein).

    Methods:
        __repr__(): Provides a formatted string representation of the FastaRecord, which includes
                    the identifier, description, and the sequence data.

    Usage example:
        record = FastaRecord('seq1', 'example description', 'AGTCAGTC')
        print(record)

    Output:
        id = seq1
        description = example description
        sequence = AGTCAGTC
    """

    def __init__(self, id, description, seq):
        self.id = id
        self.description = description
        self.seq = seq

    def __repr__(self):
        return f'id = {self.id} \n description = {self.description} \n sequence = {self.seq}'