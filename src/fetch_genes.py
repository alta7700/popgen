import argparse
import random
import time
from collections import defaultdict
from typing import TypedDict

from Bio import Entrez, SeqIO


class GroupCFG(TypedDict, total=True):
    term: str
    target_n: int
    min_hosts: int

# -----------------------------
# Поисковые группы и задания
# -----------------------------
YEARS = ["2019", "2020", "2021", "2022", "2023", "2024", "2025"]

GROUPS: dict[str, GroupCFG] = {}

def add_yearly_groups(prefix: str, base_term: str, target_n: int, min_hosts: int):
    for y in YEARS:
        name = f"{prefix}_{y}"
        term = f"{base_term} AND {y}[pdat]"
        GROUPS[name] = GroupCFG(
            term=term,
            target_n=target_n,
            min_hosts=min_hosts,
        )

# 1. SARS-CoV-2 человека
add_yearly_groups(
    "sarscov2_human",
    'SARS-CoV-2[Organism] AND human[Host] AND "complete genome"[Title]',
    target_n=4,   # можно меньше, иначе будет слишком много изолятов
    min_hosts=1,
)

# 2. Betacoronavirus (кроме SARS-CoV-2) — летучие мыши
add_yearly_groups(
    "betacoronavirus_bat",
    'Betacoronavirus[Organism] NOT SARS-CoV-2[Organism] '
    'AND bat[Host] AND "complete genome"[Title]',
    target_n=2,
    min_hosts=1,
)

# 3. Betacoronavirus (кроме SARS-CoV-2) — другие хозяева
add_yearly_groups(
    "betacoronavirus_other",
    'Betacoronavirus[Organism] NOT SARS-CoV-2[Organism] '
    'AND NOT bat[Host] AND "complete genome"[Title]',
    target_n=4,
    min_hosts=2,
)

# 4. SARS-CoV-1 человека
add_yearly_groups(
    "sarscov1_human",
    'SARS coronavirus[Organism] AND human[Host] AND "complete genome"[Title]',
    target_n=2,
    min_hosts=1,
)

# 5. MERS-CoV человека
add_yearly_groups(
    "merscov_human",
    'Middle East Respiratory Syndrome coronavirus[Organism] '
    'AND human[Host] AND "complete genome"[Title]',
    target_n=2,
    min_hosts=1,
)

# 6. Alphacoronavirus — люди
add_yearly_groups(
    "alphacoronavirus_human",
    'Alphacoronavirus[Organism] AND human[Host] AND "complete genome"[Title]',
    target_n=3,
    min_hosts=1,
)

# 7. Alphacoronavirus — прочие
add_yearly_groups(
    "alphacoronavirus_other",
    'Alphacoronavirus[Organism] AND NOT human[Host] AND "complete genome"[Title]',
    target_n=3,
    min_hosts=2,
)


# -----------------------------
# Вспомогательные функции
# -----------------------------
def esearch_ids(term: str, retmax: int = 400) -> list[str]:
    """Поиск ID записей в NCBI nuccore по esearch."""
    handle = Entrez.esearch(db="nuccore", term=term, retmax=retmax)
    result = Entrez.read(handle)
    handle.close()
    return result.get("IdList", [])


def fetch_genbank_record_by_id(nuccore_id: str):
    """Скачать одну GenBank-запись по внутреннему NCBI ID."""
    handle = Entrez.efetch(db="nuccore", id=nuccore_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record


def get_host_from_record(record) -> str:
    """Пытаемся извлечь host из source-фичи."""
    for feature in record.features:
        if feature.type == "source":
            hosts = feature.qualifiers.get("host", [])
            if hosts:
                return hosts[0]
    return "unknown"


def extract_genes_from_record(record, target_genes: set[str], group_name: str):
    """
    Вырезать указанные гены (по имени gene) из GenBank-записи.
    Возвращает список SeqRecord.
    """
    results = []
    for feature in record.features:
        if feature.type != "CDS":
            continue

        gene_names = feature.qualifiers.get("gene", [])
        if not gene_names:
            continue

        gene_name = gene_names[0]

        if gene_name not in target_genes:
            continue

        seq = feature.extract(record.seq)

        # формируем отдельный SeqRecord для FASTA
        seq_record = SeqIO.SeqRecord(
            seq,
            id=f"{record.id}|{gene_name}|{group_name}",
            description=(
                f"{gene_name} from {record.annotations.get('organism', '')}; "
                f"group={group_name}"
            ),
        )
        results.append(seq_record)

    return results


def select_records_for_group(
    group_name: str,
    config: GroupCFG,
    target_genes: set[str],
    delay: float,
) -> tuple[list[SeqIO.SeqRecord], dict[str, str]]:
    """
    Для одной группы:
    - ищем подходящие ID
    - скачиваем записи
    - выбираем нужное количество с учётом min_hosts
    - вырезаем гены
    Возвращаем:
      - список SeqRecord с генами
      - словарь {record.id: host}
    """
    print(f"\n=== Группа: {group_name} ===")
    term = config["term"]
    target_n = config["target_n"]
    max_n = target_n * 2
    min_hosts = config["min_hosts"]

    print(f"Поисковый запрос:\n  {term}")
    ids = esearch_ids(term)
    print(f"Найдено ID: {len(ids)}")
    # Перемешиваем порядок ID, чтобы не идти просто сверху вниз
    random.shuffle(ids)

    if not ids:
        print("  Ничего не найдено по этому запросу.")
        return [], {}

    selected_hosts: set[str] = set()
    host_by_record: dict[str, str] = {}
    gene_records: list[SeqIO.SeqRecord] = []

    for i, nid in enumerate(ids, start=1):
        # Останавливаемся, когда:
        # - набрали хотя бы target_n
        # - набрали не больше max_n
        # - и при этом выполняется условие для min_hosts

        if len(host_by_record) >= target_n:
            if len(host_by_record) >= max_n:
                break
            if min_hosts <= 1 or len(selected_hosts) >= min_hosts:
                # Мы уже достигли минимума и выполнили условие по хозяевам → можно остановиться
                break

        print(f"  [{i}/{len(ids)}] Скачиваем запись ID={nid} ...")
        try:
            record = fetch_genbank_record_by_id(nid)
        except Exception as e:
            print(f"    Ошибка при скачивании: {e}")
            time.sleep(delay)
            continue

        host = get_host_from_record(record)
        host_by_record[record.id] = host

        print(f"    Organism: {record.annotations.get('organism', 'unknown')}")
        print(f"    Host: {host}")

        if host != "unknown":
            selected_hosts.add(host)

        # Вырезаем гены
        records_here = extract_genes_from_record(record, target_genes, group_name)
        if not records_here:
            print(f"    Не найдено ни одного из генов {target_genes} в {record.id}")
        else:
            for gr in records_here:
                print(f"    Ген {gr.id.split('|')[1]}: длина {len(gr.seq)} нт")
            gene_records.extend(records_here)

        time.sleep(delay)

    print(
        f"Итог для группы {group_name}: {len(host_by_record)} геномов, "
        f"уникальных хозяев: {len(set(host_by_record.values()))}"
    )
    return gene_records, host_by_record


# -----------------------------
# Основной pipeline загрузки
# -----------------------------
def main():

    Entrez.email = "alta7700@mail.ru"
    Entrez.tool = "corona_auto_pipeline"

    parser = argparse.ArgumentParser(
        description=(
            "Автоматический поиск коронавирусных геномов в NCBI, "
            "выделение указанных генов и сохранение в FASTA."
        )
    )
    parser.add_argument(
        "--genes",
        nargs="+",
        help="Список имён генов, которые нужно вытащить. "
             "ORF1ab можно просто не указывать.",
    )
    parser.add_argument(
        "--output",
        help="Имя выходного FASTA-файла.",
    )
    args = parser.parse_args()
    output = args.output
    target_genes = set(args.genes)
    print(f"Будем искать гены: {', '.join(target_genes)}")

    all_gene_records: list[SeqIO.SeqRecord] = []
    hosts_summary: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))

    for group_name, cfg in GROUPS.items():
        records, host_map = select_records_for_group(
            group_name, cfg, target_genes, delay=0.5,
        )
        all_gene_records.extend(records)

        for rec_id, host in host_map.items():
            hosts_summary[group_name][host] += 1

    if not all_gene_records:
        print("Не удалось извлечь ни одной последовательности. "
              "Возможно, нужно подправить поисковые термы или имена генов.")
        return

    # Сохраняем общий FASTA
    SeqIO.write(all_gene_records, output, "fasta")
    print(f"\nГотово! Сохранено {len(all_gene_records)} последовательностей в {output}\n")

    # Краткая сводка по хозяевам
    print("Сводка по хозяевам (host) по группам:")
    for group_name, hosts in hosts_summary.items():
        print(f"  {group_name}:")
        for host, cnt in hosts.items():
            print(f"    {host}: {cnt} геном(ов)")


if __name__ == "__main__":
    main()
