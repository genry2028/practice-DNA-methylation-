import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq


class SeqDNA:
    """
    Класс для анализа ДНК.
    """

    def read_file(self, path_file):
        """
        Функция читает файл.FASTA и возвращает интерируемый объект,
        содержащий данные всех хромосом в файле.
        """
        df_list = []
        for i in SeqIO.parse(path_file, 'fasta'):
            df_list.append(self.context_identification(i.name, i.seq.upper()))
        data_frame = pl.concat(df_list, how='vertical').write_csv('DNA5.tsv')
        return data_frame

    @staticmethod
    def context_identification(name_chr, sequence) -> tuple:
        """Поиск цитозинов и контекстов метилирования ДНК."""
        position = []
        context = []
        trinuc = []
        chromosome = []
        strand = []

        for i in range(0, len(sequence) - 3):
            if sequence[i] == 'C':
                if sequence[i + 1] == 'G':
                    context.append('CG')
                else:
                    if sequence[i + 2] == 'G':
                        context.append('CHG')
                    else:
                        context.append('CHH')
                position.append(i)
                trinuc.append(str(sequence[i: i + 3]))
                chromosome.append(name_chr)
                strand.append('+')
            if sequence[i + 2] == 'G':
                if sequence[i + 1] == 'C':
                    context.append('CG')
                else:
                    if sequence[i] == 'C':
                        context.append('CHG')
                    else:
                        context.append('CHH')
                position.append(i)
                trinucletid = Seq(sequence[i: i + 3]).reverse_complement()
                trinuc.append(str(trinucletid))
                chromosome.append(name_chr)
                strand.append('-')
        return pl.DataFrame(
            {
                'chr': chromosome,
                'position': position,
                'strand': strand,
                'context': context,
                'trinuc': trinuc
            }
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
        return SeqDNA().read_file(path)


if __name__ == '__main__':
    main('exampleDNA.fa')
