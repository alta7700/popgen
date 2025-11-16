import argparse
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO


def make_short_fasta(input_path: Path) -> tuple[Path, Path, dict[str, str]]:
    """Создаёт FASTA с короткими ID и mapping-файл."""
    records = list(SeqIO.parse(str(input_path), "fasta"))
    if not records:
        print(f"Нет записей в {input_path}", file=sys.stderr)
        sys.exit(1)

    short_path = input_path.with_suffix(".short.fasta")
    map_path = input_path.with_suffix(".mapping.tsv")

    mapping: dict[str, str] = {}
    new_records = []

    for i, r in enumerate(records, start=1):
        new_id = f"seq{i}"
        mapping[new_id] = r.id
        r.id = new_id
        r.name = new_id
        r.description = ""  # чтобы gblocks не видел длинный description
        new_records.append(r)

    SeqIO.write(new_records, str(short_path), "fasta")

    with map_path.open("w") as f:
        for short_id, orig_id in mapping.items():
            f.write(f"{short_id}\t{orig_id}\n")

    print(f"[INFO] Создан сокращённый FASTA: {short_path}")
    print(f"[INFO] Mapping сохранён в: {map_path}")
    return short_path, map_path, mapping


def run_gblocks(short_fasta: Path, gblocks_cmd: str, seq_type: str):
    """Запускает gblocks на сокращённом FASTA и не паникует по коду выхода."""
    if seq_type not in ("d", "p"):
        raise ValueError("seq_type должен быть 'd' (DNA/RNA) или 'p' (protein)")

    cmd = [gblocks_cmd, str(short_fasta), f"-t={seq_type}"]
    print(f"[INFO] Запуск: {' '.join(cmd)}")

    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"[WARN] gblocks вернул код {result.returncode}, "
              f"но продолжаем, если есть trimmed-файл.")


def restore_ids(short_fasta_gb: Path, mapping: dict[str, str], output_path: Path):
    """Читает trimmed FASTA от gblocks и восстанавливает исходные ID."""
    records = list(SeqIO.parse(str(short_fasta_gb), "fasta"))
    if not records:
        print(f"[ERROR] Нет записей в {short_fasta_gb}", file=sys.stderr)
        sys.exit(1)

    restored = []
    for r in records:
        short_id = r.id
        if short_id not in mapping:
            print(f"[WARN] ID {short_id} не найден в mapping, оставляю как есть")
            restored.append(r)
            continue

        orig_id = mapping[short_id]
        r.id = orig_id
        r.name = orig_id
        r.description = ""
        restored.append(r)

    SeqIO.write(restored, str(output_path), "fasta")
    print(f"[INFO] Итоговый файл с восстановленными ID: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Обёртка для gblocks: сокращает имена, запускает gblocks, "
            "а затем восстанавливает исходные ID."
        )
    )
    parser.add_argument(
        "input",
        help="Входной FASTA-файл с выравниванием",
    )
    parser.add_argument(
        "--gblocks-cmd",
        default="gblocks",
        help="Имя/путь к команде gblocks (по умолчанию: gblocks)",
    )
    parser.add_argument(
        "--seq-type",
        default="d",
        choices=["d", "p"],
        help="Тип последовательностей: d = DNA/RNA, p = protein (по умолчанию: d)",
    )
    parser.add_argument(
        "--keep-intermediate",
        action="store_true",
        help="Не удалять *.short.fasta и файл от gblocks (по умолчанию удаляются).",
    )

    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"[ERROR] Файл не найден: {input_path}", file=sys.stderr)
        sys.exit(1)

    # 1. Делаем короткие имена и mapping
    short_fasta, map_path, mapping = make_short_fasta(input_path)

    # 2. Запускаем gblocks
    run_gblocks(short_fasta, args.gblocks_cmd, args.seq_type)

    # 3. Gblocks создаёт файл с суффиксом '-gb'
    short_gb = Path(str(short_fasta) + "-gb")
    if not short_gb.exists():
        print(f"[ERROR] Ожидался файл {short_gb}, но его нет. "
              f"Проверь вывод gblocks.", file=sys.stderr)
        sys.exit(1)

    # 4. Восстанавливаем ID
    # имя вида: corona_S_aligned.gblocks.fasta
    tmp = input_path.with_suffix("")          # corona_S_aligned
    output_path = tmp.with_suffix(".gblocks.fasta")

    restore_ids(short_gb, mapping, output_path)

    # 5. При желании удаляем промежуточные файлы
    if not args.keep_intermediate:
        try:
            short_fasta.unlink()
            map_path.unlink()
            # html отчёт gblocks (если есть) оставляем
            # сам trimmed файл тоже можно удалить, он уже восстановлен
            short_gb.unlink()
        except OSError:
            pass

    print("[INFO] Готово.")


if __name__ == "__main__":
    main()
