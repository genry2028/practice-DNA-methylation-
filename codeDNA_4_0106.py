import polars as pl
from Bio import SeqIO


class SeqDNA:
    """
    Класс для анализа ДНК.
    """

    __slots__ = ['path_file', 'position', 'context', 'trinuc', 'chromosome']

    def __init__(self, path_file) -> str:
        self.path_file = path_file
        self.position = []
        self.context = []
        self.trinuc = []
        self.chromosome = []

    # Не совсем правильно, что функция, с данными которой пользователь работать  непосредственно (а не приватная для
    # класса) возвращает свои же аттрибуты. В случае списков, словарей, классов и т.д. Python возвращает указатель на
    # объект, и изменение одного объекта повлечет и изменение исходного. Ниже – пример:
    # class Example:
    #     def __init__(self):
    #         self.var1 = [1, 2, 3]
    #
    #     def return_var1(self):
    #         return self.var1
    #
    # example = Example()
    # var2 = example.return_var1()
    #
    # example.var1
    # >>> [1, 2, 3]
    # var2
    # >>> [1, 2, 3]
    # var2.pop()
    # >>> var2
    # [1, 2]
    # example.var1
    # >>> [1, 2]
    # Как видим, обе переменных поменялись, хотя меняли только var2.
    #
    # Если нам все-таки нужно возвратить, как в данном случае, свои аттрибуты, но без возможности их изменения, то на
    # этот случай есть стандартная библиотека copy.
    #
    # САМ КОММЕНТАРИЙ :) – Ну а конкретно здесь нам гораздо проще было бы возвращать сразу DataFrame.
    def read_file(self):
        """
        Функция читает файл.FASTA и возвращает интерируемый объект,
        содержащий данные всех хромосом в файле.
        """
        records = [record for record in SeqIO.parse(self.path_file, 'fasta')]
        for i in records:
            self.get_analisis(i.name, i.seq.upper())
        return (
            self.position, self.context, self.trinuc, self.chromosome
        )

    # Попробуйте переименовать, чтобы было лучше понятно предназначение функции
    # (если пользуйтесь PyCharm – Shift+F6 - очень удобно!)
    #
    # Здесь тоже можно возвращать DataFrame. Посмотрите ф-ию polars.concat()
    def get_analisis(self, name_chr, sequence) -> tuple:
        """Анализ последовательности."""
        # Посмотрите на последовательности, которые соответствуют контексты. Есть ли в них что-то общее, что позволит
        # оптимизировать определение контекста? Потому что на данный момент требуется произвести в худшем случае
        # 16 сравнений. (Можно обойтись простыми if...elif...else)
        #
        # P.S. Определять контекст с помощью словаря – хорошая идея, однако если идти данным методом, то лучше было бы
        # сначала составить DataFrame только с координатами и тринуклеотидом, а затем создать столбец контекстов путем
        # встроенного функционала polars.
        # full_df = trinuc_df.with_columns(pl.col("trinuc").replace(trinuc_to_context_dict).alias("context"))
        # Здесь trinuc_to_context_dict - словарь по типу {'CGA': 'CG', 'CAA': 'CHH'...}
        options_H = {
            'CG': ('CGA', 'CGC', 'CGT', 'CGG'),
            'CHG': ('CAG', 'CCG', 'CTG'),
            'CHH': (
                'CAA', 'CAT', 'CAC', 'CTA', 'CTT', 'CTC', 'CCA', 'CCT', 'CCC'
                )
        }
        for i in range(0, len(sequence) - 3):
            if sequence[i] == 'C':
                site = sequence[i: i + 3]
                for key, value in options_H.items():
                    if site in value:
                        self.position.append(i)
                        self.context.append(key)
                        self.trinuc.append(str(site))
                        self.chromosome.append(name_chr)
        return (
            self.position, self.context, self.trinuc, self.chromosome
            )


# Отрицательная нить ДНК здесь никак не учитывается.
def main(path=None, strand_DNA='+'):
    """Основная функция, возвращающая таблицу данных."""
    if path is None:
        print(
            'Для запуска программы и начала анализа '
            'введите адрес файла в формате Fasta.'
        )
        path = input()
        main(path)
    else:
        position, context, trinuc, chromosome = SeqDNA(path).read_file()
        # Вот это убрать в read_file
        data_frame = pl.DataFrame(
            {
                'chr': chromosome,
                'position': position,
                'strand': [strand_DNA for k in range(len(trinuc))],
                'context': context,
                'trinuc': trinuc
            }
        )
        return data_frame.write_csv('for3_arabi_ict_DNA.tsv')


if __name__ == '__main__':
    main()