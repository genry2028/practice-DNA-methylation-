import itertools

import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq


class SeqDNA:
    """
    Класс для анализа ДНК.
    """
    # Всегда когда данные предсказуемы лучше явно указывать схему DataFrame, иначе polars пытается угадать ее из входных
    # данных, что может повлечь за собой неожиданные несоответствия типов двух, казалось бы, похожих таблиц. К тому же
    # это повышает быстродействие.
    # Здесь мы указываем схему как аттрибут класса, не требующий его инициализации.
    trinuc_schema = dict(
        position = pl.UInt64,
        context=pl.String,
        trinuc=pl.String,
        chromosome=pl.String,
        strand=pl.Boolean
    )

    def read_file(self, path_file):
        """
        Функция читает файл.FASTA и возвращает интерируемый объект,
        содержащий данные всех хромосом в файле.
        """
        # Здесь мы создаем пустую DataFrame с указанной схемой.
        df = pl.DataFrame(schema=self.trinuc_schema)

        for i in SeqIO.parse(path_file, 'fasta'):
            # Здесь чисто для того, чтобы убрать пару строчек изменение.
            # В polars есть три метода:
            # pl.concat – объединяет список таблиц в одну. 
            # pl.DataFrame.vstack – не изменяет исходную таблицу, но возвращает расширенную
            # pl.DataFrame.extend – расширяет исходную таблицу, наподобие list.append
            df.extend(self.chrom_trinuc_df(i.name, i.seq.upper()))

        return df
    
    # Сначала см. комментарий ниже.
    # Это просто обертка над extract_trinuc, т.к. numba не сможет DataFrame скомпилировать, этот пакет умеет работать со
    # стандартными библиотеками и numpy, но тем не менее порой позволяет ускорить работу в 2-3 раза.
    def chrom_trinuc_df(self, name_chr, sequence):
        return (
            pl.DataFrame(
                self.extract_trinuc(sequence) + tuple(itertools.repeat(name_chr, len(sequence))),
                schema=self.trinuc_schema
            )
        )

    # Отличное решение сделать эту функцию @staticmethod, т.к. теперь она независима от класса и ее можно, например,
    # скомпилировать (см. пакет numba) для повышения быстродействия.
    # P.S. Прошу прощения, что переименовал... identification – длинное слово :(
    @staticmethod
    def extract_trinuc(sequence) -> tuple:
        """Поиск цитозинов и контекстов метилирования ДНК."""
        position = []
        context = []
        trinuc = []
        # Убрал здесь создание списка хромосом, поскольку он и так одинаковый, и мало ли, потребуется просто
        # проанализировать последовательность без учета хромосомы.
        strand = []

        for i in range(0, len(sequence) - 3):
            if sequence[i] == 'C':
                if sequence[i + 1] == 'G':
                    context.append('CG')
                elif sequence[i + 2] == 'G':
                    context.append('CHG')
                else:
                    context.append('CHH')
                position.append(i)
                trinuc.append(str(sequence[i: i + 3]))
                strand.append(True)
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
                strand.append(False)

        return position, context, trinuc, strand


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
