import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq


# Сейчас никаких грубых ошибок уже нет, поэтому комментарии будут касаться, можно сказать, принципов программирования.
#
# Исходная задача стояла в написании класса, который неким образом работает с геномной последовательностью, в нашем
# случае – извлекает контекст метилирования. Но ведь в теории, например, мы могли бы расширить данный класс, чтобы он
# считал количество нуклеотидов или другую статистику на основании последовательности. Вопрос: насколько Ваш класс
# пригоден для масштабирования?
#
# Ниже будут комментарии, которые ответят на вышепоставленный вопрос.
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
        # Все эти списки лучше не создавать здесь.
        self.position_plus = []
        self.position_minus = []
        self.context_plus = []
        self.context_minus = []
        self.trinuc_plus = []
        self.trinuc_minus = []
        self.chromosome_plus = []
        self.chromosome_minus = []

    # Вот функция называется read_file, а файла в качестве параметра не получает. Возникает дилема – либо класс называть
    # SeqFile, а метод read, и позиционировать класс – как класс для работы с *файлом* последовательности, либо убирать
    # path_file из аттрибутов и добавлять его в параметр этой функции (лично я бы сделал 2 вариант).
    def read_file(self):
        """
        Функция читает файл.FASTA и возвращает интерируемый объект,
        содержащий данные всех хромосом в файле.
        """
        # Зачем хранить в памяти лишние последовательности и сразу читать все хромосомы, если мы можем итерировать сразу?
        # records = [record for record in SeqIO.parse(self.path_file, 'fasta')]
        # for i in records:

        for i in SeqIO.parse(self.path_file, 'fasta'):
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
        # Пусть лучше context_identification возвращает таблицу – сразу и +, и - нити в ней будут. Потом список из
        # таблиц для каждой хромосомы соединим методом pl.concat

        df = pl.concat(
            [data_frame1, data_frame2], how='vertical'
        ).write_csv('DNA_5.tsv')
        return df

    # Чем меньше функция имеет доступа к внешним для себя объектам – тем лучше! Данную функцию в принципе можно было бы
    # сделать статичной (@staticmethod), ведь если бы все эти списки, которые создаются у Вас в __init__, создавались бы
    # здесь, то функции бы вообще не потребовалось знать, что там происходит в остальном классе, а это всегда упрощает
    # отладку.
    def context_identification(self, name_chr, sequence) -> tuple:
        """Поиск цитозинов и контекстов метилирования ДНК."""
        for i in range(0, len(sequence) - 3):
            if sequence[i] == 'C':
                # Посмотрите на условия. Можно ли их расположить в другом порядке, чтобы сравнений было меньше?
                if sequence[i + 1] != 'G' and sequence[i + 2] != 'G':
                    self.context_plus.append('CHH')
                elif sequence[i + 1] != 'G' and sequence[i + 2] == 'G':
                    self.context_plus.append('CHG')
                else:
                    self.context_plus.append('CG')
                self.position_plus.append(i)
                self.trinuc_plus.append(str(sequence[i: i + 3]))
                self.chromosome_plus.append(name_chr)
            # Отлично!
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
        # Почему бы здесь нам не возвращать DataFrame сразу?
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
