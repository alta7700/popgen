import argparse
from pathlib import Path

from Bio import SeqIO


def guess_gene_name(record, allowed_genes):
    """
    Пытаемся определить имя гена из ID/description.
    Ориентируемся в первую очередь на формат ACCESSION|GENE|GROUP.
    """
    # 1. Основной случай: ACCESSION|GENE|GROUP
    parts = record.id.split("|")
    if len(parts) >= 2:
        gene = parts[1]
        if gene in allowed_genes:
            return gene

    # 2. На всякий случай попробуем поискать имя гена в description
    for g in allowed_genes:
        if f"|{g}|" in record.description or f" {g} " in record.description:
            return g

    return None


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Разделение общего FASTA-файла на отдельные файлы по генам "
            "(например, S и N)."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Входной FASTA-файл (результат предыдущего пайплайна).",
    )
    parser.add_argument(
        "--genes",
        nargs="+",
        required=True,
        help="Имена генов, которые нужно выделить.",
    )
    parser.add_argument(
        "--prefix",
        help="Префикс для имён выходных файлов.",
    )

    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Входной файл не найден: {input_path}")

    allowed_genes = set(args.genes)
    records_by_gene = {g: [] for g in allowed_genes}

    print(f"Читаем {input_path} ...")
    total = 0
    skipped = 0

    for record in SeqIO.parse(str(input_path), "fasta"):
        total += 1
        gene = guess_gene_name(record, allowed_genes)
        if gene is None:
            skipped += 1
            print(f"  [!] Не удалось определить ген для записи: {record.id}")
            continue

        records_by_gene[gene].append(record)

    print(f"\nВсего записей: {total}")
    print(f"Пропущено (неопознанный ген): {skipped}")

    # Пишем по файлам
    for gene, recs in records_by_gene.items():
        if not recs:
            print(f"  [!] Для гена {gene} не найдено ни одной последовательности.")
            continue
        out_path = input_path.with_name(f"{args.prefix}_{gene}.fasta")
        SeqIO.write(recs, out_path, "fasta")
        print(f"  Ген {gene}: {len(recs)} последовательностей → {out_path}")

    print("\nГотово.")


if __name__ == "__main__":
    main()
