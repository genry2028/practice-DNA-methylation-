import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq


class SeqDNA:
    """
    Класс для анализа ДНК.
    """

    __slots__ = [
        'path_file', 'position_plus', 'context_plus', 'trinuc_plus',
        'chromosome_plus', 'position_minus', 'context_minus',
        'trinuc_minus', 'chromosome_minus'
    ]

    def __init__(self, path_file) -> str:
        self.path_file = path_file
        self.position_plus = []
        self.position_minus = []
        self.context_plus = []
        self.context_minus = []
        self.trinuc_plus = []
        self.trinuc_minus = []
        self.chromosome_plus = []
        self.chromosome_minus = []

    def read_file(self):
        """
        Функция читает файл.FASTA и возвращает интерируемый объект,
        содержащий данные всех хромосом в файле.
        """
        records = [record for record in SeqIO.parse(self.path_file, 'fasta')]
        for i in records:
            self.context_identification(i.name, i.seq.upper())
        data_frame1 = pl.DataFrame(
            {
                'chr': self.chromosome_plus,
                'position': self.position_plus,
                'strand': ['+' for k in range(len(self.trinuc_plus))],
                'context': self.context_plus,
                'trinuc': self.trinuc_plus
            }
        )
        data_frame2 = pl.DataFrame(
            {
                'chr': self.chromosome_minus,
                'position': self.position_minus,
                'strand': ['-' for k in range(len(self.trinuc_minus))],
                'context': self.context_minus,
                'trinuc': self.trinuc_minus
            }
        )
        df = pl.concat(
            [data_frame1, data_frame2], how='vertical'
        ).write_csv('DNA_5.tsv')
        return df

    def context_identification(self, name_chr, sequence) -> tuple:
        """Поиск цитозинов и контекстов метилирования ДНК."""
        for i in range(0, len(sequence) - 3):
            if sequence[i] == 'C':
                if sequence[i + 1] != 'G' and sequence[i + 2] != 'G':
                    self.context_plus.append('CHH')
                elif sequence[i + 1] != 'G' and sequence[i + 2] == 'G':
                    self.context_plus.append('CHG')
                else:
                    self.context_plus.append('CG')
                self.position_plus.append(i)
                self.trinuc_plus.append(str(sequence[i: i + 3]))
                self.chromosome_plus.append(name_chr)
            if sequence[i + 2] == 'G':
                if sequence[i] != 'C' and sequence[i + 1] != 'C':
                    self.context_minus.append('CHH')
                elif sequence[i + 1] != 'C' and sequence[i] == 'C':
                    self.context_minus.append('CHG')
                else:
                    self.context_minus.append('CG')
                self.position_minus.append(i)
                trinuc = Seq(sequence[i: i + 3]).reverse_complement()
                self.trinuc_minus.append(str(trinuc))
                self.chromosome_minus.append(name_chr)
        return (
            self.position_plus.copy(), self.context_plus.copy(),
            self.trinuc_plus.copy(), self.chromosome_plus.copy(),
            self.position_minus.copy(), self.context_minus.copy(),
            self.trinuc_minus.copy(), self.chromosome_minus.copy()
        )


def main(path=None):
    """Основная функция, возвращающая таблицу данных."""
    if path is None:
        print(
            'Для запуска программы и начала анализа '
            'введите адрес файла в формате Fasta.'
        )
        path = input()
        main(path)
    else:
        return SeqDNA(path).read_file()


if __name__ == '__main__':
    main()
