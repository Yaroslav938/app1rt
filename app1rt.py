import streamlit as st
import re
import pandas as pd
import io
import plotly.express as px
import plotly.graph_objects as go
import math
from io import BytesIO

try:
    from striprtf.striprtf import rtf_to_text
except ImportError:
    rtf_to_text = None

# ─── Настройка страницы ───
st.set_page_config(
    page_title="RT-PCR Анализатор",
    page_icon="🧬",
    layout="wide",
)

st.markdown("""
<style>
    .main .block-container { max-width: 1200px; padding-top: 2rem; }
    .stTabs [data-baseweb="tab-list"] { gap: 8px; }
    .stTabs [data-baseweb="tab"] { padding: 10px 24px; font-weight: 600; }
    div[data-testid="stMetricValue"] { font-size: 1.1rem; }
    .result-positive { color: #d32f2f; font-weight: bold; }
    .result-negative { color: #388e3c; font-weight: bold; }
    .result-warning  { color: #f57c00; font-weight: bold; }
    .result-info     { color: #1976d2; font-weight: bold; }
    .header-card {
        background: linear-gradient(135deg, #1a237e 0%, #283593 100%);
        color: white; padding: 1.5rem 2rem; border-radius: 12px;
        margin-bottom: 1.5rem;
    }
    .header-card h1 { color: white; margin: 0 0 0.3rem 0; font-size: 1.8rem; }
    .header-card p  { color: #bbdefb; margin: 0; font-size: 0.95rem; }
    .info-box { background: #e3f2fd; border-left: 4px solid #1976d2; padding: 12px 16px; border-radius: 0 8px 8px 0; margin: 8px 0; }
    .warn-box { background: #fff3e0; border-left: 4px solid #f57c00; padding: 12px 16px; border-radius: 0 8px 8px 0; margin: 8px 0; }
    .error-box { background: #ffebee; border-left: 4px solid #d32f2f; padding: 12px 16px; border-radius: 0 8px 8px 0; margin: 8px 0; }
    .success-box { background: #e8f5e9; border-left: 4px solid #388e3c; padding: 12px 16px; border-radius: 0 8px 8px 0; margin: 8px 0; }
</style>
""", unsafe_allow_html=True)


# ═══════════════════════════════════════════════════════════════════
# БЛОК 1: УТИЛИТЫ ДЛЯ ПАРСИНГА И ЭКСПОРТА ФАЙЛОВ
# ═══════════════════════════════════════════════════════════════════

def _try_decode(content_bytes: bytes) -> str:
    """Умная декодировка с приоритетом cp1251 (РФ-протоколы)."""
    for enc in ("cp1251", "utf-8", "latin-1"):
        try:
            text = content_bytes.decode(enc)
            if any(kw in text for kw in ("Номер", "Дата", "Оператор", "лунки", "Fam", "Hex", "Rox")):
                return text
        except UnicodeDecodeError:
            continue
    return content_bytes.decode("cp1251", errors="ignore")


def parse_rtf_protocol(content_bytes: bytes) -> dict:
    """
    Парсит RTF-отчёт амплификатора ДНК-Технологии.
    """
    raw = _try_decode(content_bytes)

    if rtf_to_text is not None and raw.strip().startswith("{\\rtf"):
        text = rtf_to_text(raw)
    else:
        text = raw

    result = {
        "date": "Неизвестно",
        "operator": "Неизвестно",
        "protocol_number": "Неизвестно",
        "result_file": "Неизвестно",
        "test_name": "Неизвестно",
        "amplification_program": "Неизвестно",
        "samples": [],
    }

    lines = text.splitlines()
    sample_section = False

    ch_indices = {}   
    well_idx = -1     
    name_idx = -1     

    for line in lines:
        lc = line.strip()
        if not lc:
            continue

        parts_meta = [p.strip() for p in lc.split("|")]

        def _extract_meta(keyword):
            for i, p in enumerate(parts_meta):
                if keyword in p:
                    for j in range(i + 1, len(parts_meta)):
                        if parts_meta[j]:
                            return parts_meta[j]
            return None

        if "Дата:" in lc:
            v = _extract_meta("Дата:")
            if v: result["date"] = v
        if "Оператор:" in lc:
            v = _extract_meta("Оператор:")
            if v: result["operator"] = v
        if "Номер протокола:" in lc:
            v = _extract_meta("Номер протокола:")
            if v: result["protocol_number"] = v
        if "Файл с результатами:" in lc or "Файл:" in lc:
            v = _extract_meta("Файл")
            if v and "Файл" not in v: result["result_file"] = v
        if "Тест:" in lc:
            v = _extract_meta("Тест:")
            if v: result["test_name"] = v
        if "Программа амплификации:" in lc:
            v = _extract_meta("Программа амплификации:")
            if v: result["amplification_program"] = v

        lc_lower = lc.lower()
        if ("номер лунки" in lc_lower or "номер_лунки" in lc_lower) and (
                "идентификатор" in lc_lower or "пробирка" in lc_lower or "образец" in lc_lower):
            sample_section = True
            parts = [p.strip() for p in lc.split("|")]
            for i, p in enumerate(parts):
                pl = p.lower()
                if "номер лунки" in pl or "номер_лунки" in pl:
                    well_idx = i
                elif "идентификатор" in pl or "пробирка" in pl or "образец" in pl:
                    name_idx = i
                elif "fam" in pl or "green" in pl:
                    ch_indices["Fam"] = i
                elif "hex" in pl or "yellow" in pl:
                    ch_indices["Hex"] = i
                elif "rox" in pl or "orange" in pl:
                    ch_indices["Rox"] = i
                elif "cy5.5" in pl:
                    ch_indices["Cy5.5"] = i
                elif "cy5" in pl or "red" in pl:
                    ch_indices["Cy5"] = i
            continue

        if sample_section:
            if "* " in lc or ("|*" in lc):
                sample_section = False
                continue

            parts = [p.strip() for p in lc.split("|")]

            well = ""
            found_well_data_i = -1
            for i, p in enumerate(parts):
                if re.match(r'^[A-H]\d{1,2}$', p):
                    well = p
                    found_well_data_i = i
                    break

            if not well:
                continue

            offset = found_well_data_i - well_idx if well_idx >= 0 else 0

            name = ""
            name_data_i = name_idx + offset if name_idx >= 0 else found_well_data_i + 1
            if 0 <= name_data_i < len(parts) and parts[name_data_i]:
                name = parts[name_data_i]
            elif found_well_data_i + 1 < len(parts):
                name = parts[found_well_data_i + 1]
            name = re.sub(r'\s*\(.*?\)\s*$', '', name).strip()

            ct_values = {}
            for ch, ci in ch_indices.items():
                actual_ci = ci + offset
                if 0 <= actual_ci < len(parts):
                    val = parts[actual_ci].replace(",", ".").strip()
                    if re.match(r'^\d+\.?\d*$', val) and float(val) > 0:
                        ct_values[ch] = float(val)
                    else:
                        ct_values[ch] = None
                else:
                    ct_values[ch] = None

            result["samples"].append({
                "well": well,
                "name": name,
                "ct": ct_values,
            })

    return result


def parse_r48(content_bytes: bytes) -> dict:
    """Парсинг сырых файлов .r48/.r96 от ДНК-Технологии."""
    text = None
    for enc in ("cp1251", "utf-8", "latin-1"):
        try:
            text = content_bytes.decode(enc)
            break
        except UnicodeDecodeError:
            pass
    if not text:
        raise ValueError("Ошибка кодировки файла данных")

    lines = text.split('\n')

    num_cycles = 0
    for line in lines[:100]:
        if 'XCYC' in line:
            m = re.search(r'XCYC\s+(\d+)', line)
            if m:
                c = int(m.group(1))
                if c > 20:
                    num_cycles = c
    if num_cycles == 0:
        num_cycles = 50

    channels = ['Fam', 'Hex', 'Rox', 'Cy5', 'Cy5.5']
    raw_data = {ch: [] for ch in channels}

    for line in lines:
        parts = line.strip().split()
        if len(parts) < 10:
            continue
        ch_raw = parts[0].split('_')[0]
        if ch_raw in channels:
            try:
                vals = [float(x) for x in parts[7:]]
                raw_data[ch_raw].append(vals)
            except ValueError:
                continue

    result = {}
    for ch in channels:
        if len(raw_data[ch]) > 0:
            if len(raw_data[ch]) >= num_cycles * 2:
                reads_per_cycle = len(raw_data[ch]) // num_cycles
                blocks = [raw_data[ch][i::reads_per_cycle][:num_cycles] for i in range(reads_per_cycle)]
                block = max(blocks, key=lambda b: sum(b[-1]))
            else:
                block = raw_data[ch][:num_cycles]

            num_wells = len(block[0])

            if num_wells == 48:
                cols, rws = 8, 6
            elif num_wells == 96:
                cols, rws = 12, 8
            else:
                cols, rws = num_wells, 1

            well_names = []
            for r in range(rws):
                for c in range(cols):
                    well_names.append(f"{chr(65 + r)}{c + 1}")

            if len(well_names) != num_wells:
                well_names = [str(i + 1) for i in range(num_wells)]

            df = pd.DataFrame(block, columns=well_names)
            df.index = range(1, len(df) + 1)
            df.index.name = 'Цикл'
            result[ch] = df

    return result


def calculate_raw_summary(raw_dfs: dict, protocol_data: dict = None) -> pd.DataFrame:
    """Создаёт итоговую таблицу Ct для вкладки сырых данных."""
    if not raw_dfs:
        return pd.DataFrame()
    first_ch = list(raw_dfs.keys())[0]
    wells = raw_dfs[first_ch].columns.tolist()

    rows = []
    display_channels = ["Fam", "Hex", "Rox", "Cy5", "Cy5.5"]

    rtf_map = {}
    if protocol_data and "samples" in protocol_data:
        for s in protocol_data["samples"]:
            rtf_map[s["well"]] = s

    for w in wells:
        sample_name = rtf_map[w]["name"] if w in rtf_map else ""
        row = {"Лунка": w, "Образец": sample_name}
        has_signal = False

        for ch in display_channels:
            ct_val = None

            if protocol_data and w in rtf_map:
                val = rtf_map[w].get("ct", {}).get(ch)
                if val is not None:
                    ct_val = float(val)
                    has_signal = True

            elif ch in raw_dfs and w in raw_dfs[ch].columns:
                series = raw_dfs[ch][w]
                if len(series) > 15:
                    smooth = series.rolling(window=3, min_periods=1, center=True).mean()
                    base_region = smooth.loc[5:15]
                    baseline = base_region.mean()
                    noise_amp = base_region.max() - base_region.min()
                    norm = smooth - baseline
                    max_peak = norm.max()
                    peak_idx = norm.idxmax()

                    if max_peak >= 10.0 and max_peak > (noise_amp * 3.0) and peak_idx > 15:
                        max_slope = norm.diff(periods=2).max()
                        final_val = norm.iloc[-3:].mean()
                        if max_slope >= 1.5 and final_val >= 5.0:
                            threshold = max(max_peak * 0.1, 5.0)
                            cross = norm[(norm.index <= peak_idx) & (norm >= threshold)]
                            if not cross.empty and cross.index[0] > 5:
                                ct_val = float(round(cross.index[0], 1))
                                has_signal = True

            row[f"{ch} Ct"] = f"{ct_val:.1f}" if ct_val is not None else "—"

        if sample_name or has_signal:
            if not row["Образец"]:
                row["Образец"] = "Неизвестный образец"
            rows.append(row)

    return pd.DataFrame(rows)


def to_excel(df: pd.DataFrame) -> bytes:
    """Экспорт DataFrame в Excel."""
    output = BytesIO()
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        df.to_excel(writer, index=False, sheet_name='Результаты')
        ws = writer.sheets['Результаты']
        for i, col in enumerate(df.columns):
            column_len = max(df[col].astype(str).map(len).max(), len(col)) + 2
            ws.set_column(i, i, column_len)
    return output.getvalue()


# ═══════════════════════════════════════════════════════════════════
# БЛОК 2: АНАЛИЗАТОРЫ — МАТЕМАТИКА ТОЧНО ПО EXCEL
# ═══════════════════════════════════════════════════════════════════

def _ct(sample: dict, ch: str):
    """Безопасное извлечение Ct (float или None)."""
    if not sample:
        return None
    return sample.get("ct", {}).get(ch)


def analyze_plant_gmo(protocol: dict, cfg: dict, is_swab: bool = False) -> dict:
    samples = protocol["samples"]
    pko_list, oko_list, test_samples = [], [], []

    for s in samples:
        n = s["name"].lower().strip()
        if "пко" in n or n.startswith("пко") or "pko" in n:
            pko_list.append(s)
        elif "око" in n or n.startswith("око") or "nko" in n:
            oko_list.append(s)
        else:
            test_samples.append(s)

    if not pko_list or not oko_list:
        return {"error": "Не найден ПКО или ОКО. Проверьте названия лунок в протоколе (должны содержать 'пко' и 'око')."}

    pko = pko_list[0]
    oko = oko_list[0]

    pko_plant  = _ct(pko, "Hex")   
    pko_35s    = _ct(pko, "Rox")   
    pko_nos    = _ct(pko, "Fam")   
    pko_vpk    = _ct(pko, "Cy5")   

    oko_plant  = _ct(oko, "Hex")
    oko_35s    = _ct(oko, "Rox")
    oko_nos    = _ct(oko, "Fam")
    oko_vpk    = _ct(oko, "Cy5")

    pko_invalid = (
        (pko_plant is None or pko_plant > cfg['pko_plant_max']) or
        (pko_35s   is None or pko_35s   > cfg['pko_prom_max']) or
        (pko_nos   is None or pko_nos   > cfg['pko_prom_max'])
    )
    pko_status = "Не принимается ПКО" if pko_invalid else "Принимается"
    pko_action = "Повторить анализ" if pko_invalid else ""

    oko_vpk_ok = (oko_vpk is not None and oko_vpk <= cfg['vpk_max'])

    oko_plant_det = _detect_plant(oko_plant, pko_plant, cfg)
    oko_35s_det   = "-"
    oko_nos_det   = "-"
    if oko_plant_det in ("+", "мало") and pko_plant is not None and pko_35s is not None and pko_nos is not None:
        oko_35s_det = _detect_promoter(oko_35s, oko_plant, pko_35s, pko_plant, oko_plant_det, pko_35s, cfg)
        oko_nos_det = _detect_promoter(oko_nos, oko_plant, pko_nos, pko_plant, oko_plant_det, pko_nos, cfg)

    oko_contaminated = (oko_plant_det == "+" or oko_35s_det == "+" or oko_nos_det == "+")

    if oko_contaminated:
        oko_status = "Контаминация"
        oko_action = "Провести деконтаминацию"
    elif not oko_vpk_ok:
        oko_status = "Не принимается ВПК"
        oko_action = "Повторить анализ"
    else:
        oko_status = "Принимается"
        oko_action = ""

    sample_results = []
    for s in test_samples:
        s_plant = _ct(s, "Hex")
        s_35s   = _ct(s, "Rox")
        s_nos   = _ct(s, "Fam")
        s_vpk   = _ct(s, "Cy5")

        if is_swab:
            plant_det = _swab_det(s_plant)
            prom_det  = _swab_det(s_35s)
            nos_det   = _swab_det(s_nos)
        else:
            plant_det = _detect_plant(s_plant, pko_plant, cfg)
            prom_det  = _detect_promoter(s_35s, s_plant, pko_35s, pko_plant, plant_det, pko_35s, cfg)
            nos_det   = _detect_promoter(s_nos, s_plant, pko_nos, pko_plant, plant_det, pko_nos, cfg)

        if s_vpk is not None:
            if oko_vpk is not None:
                vpk_det = "+" if s_vpk <= oko_vpk + cfg['vpk_delta'] else "-"
            else:
                vpk_det = "+" if s_vpk <= cfg['vpk_max'] else "-"
        else:
            vpk_det = "-"

        conclusion, action = _get_conclusion_gmo(plant_det, prom_det, nos_det, vpk_det, is_swab)

        sample_results.append({
            "name": s["name"], "well": s["well"],
            "plant_ct": s_plant, "plant_det": plant_det,
            "promoter_ct": s_35s, "promoter_det": prom_det,
            "nos_ct": s_nos, "nos_det": nos_det,
            "vpk_ct": s_vpk, "vpk_det": vpk_det,
            "conclusion": conclusion, "action": action,
        })

    return {
        "pko": {
            "name": pko["name"], "well": pko["well"],
            "plant_ct": pko_plant, "promoter_ct": pko_35s,
            "nos_ct": pko_nos, "vpk_ct": pko_vpk,
            "status": pko_status, "action": pko_action,
        },
        "oko": {
            "name": oko["name"], "well": oko["well"],
            "plant_ct": oko_plant, "vpk_ct": oko_vpk,
            "status": oko_status, "action": oko_action,
        },
        "samples": sample_results,
        "pko_plant_avg": pko_plant,
        "pko_35s_avg":   pko_35s,
        "pko_nos_avg":   pko_nos,
    }


def _detect_plant(ct_plant, pko_plant, cfg) -> str:
    if ct_plant is None or pko_plant is None:
        return "-"
    if ct_plant > pko_plant + cfg['delta_low']:
        return "-"
    elif ct_plant > pko_plant + cfg['delta_pos']:
        return "мало"
    else:
        return "+"


def _detect_promoter(ct_prom, ct_plant_sample, pko_prom, pko_plant, plant_det, pko_prom_ref, cfg) -> str:
    if ct_prom is None:
        return "-"
    if plant_det == "+":
        if ct_plant_sample is None or pko_prom is None or pko_plant is None:
            return "-"
        delta_sample = ct_prom - ct_plant_sample
        delta_pko    = pko_prom - pko_plant
        if delta_sample > delta_pko + cfg['delta_pos']:
            return "-"
        else:
            return "+"
    elif plant_det == "мало":
        if pko_prom_ref is None:
            return "-"
        if ct_prom > pko_prom_ref + cfg['delta_pos']:
            return "-"
        else:
            return "+"
    else:
        return "-"


def _swab_det(ct_val) -> str:
    if ct_val is None:
        return "-"
    if ct_val <= 30:
        return "+"
    return "мало"


def _get_conclusion_gmo(plant_det, prom_det, nos_det, vpk_det, is_swab) -> tuple:
    if is_swab:
        if plant_det == "+" or prom_det == "+" or nos_det == "+":
            return "ОБНАРУЖЕНА КОНТАМИНАЦИЯ!", "Провести повторную деконтаминацию"
        elif plant_det == "мало" or prom_det == "мало" or nos_det == "мало":
            return "Следовая контаминация (мало ДНК)", "Обратить внимание / повторная деконтаминация"
        elif vpk_det == "+":
            return "Смыв чистый (контаминации нет)", "Не требуется"
        else:
            return "Реакция невалидна (ВПК не прошёл)", "Повторить анализ"
    else:
        if plant_det == "+" and (prom_det == "+" or nos_det == "+"):
            return (
                "ГМ растение; количество ДНК выше предела количественного определения",
                "Идентификация и определение количества %",
            )
        elif plant_det == "+" and prom_det == "-" and nos_det == "-":
            return (
                "Промоторы 35S, FMV и терминатор NOS не обнаружены",
                "Обнаружение сои линий BPS-CV-127-9 и MON87701",
            )
        elif plant_det == "мало" and prom_det == "-" and nos_det == "-":
            return (
                "Промоторы 35S, FMV и терминатор NOS не обнаружены",
                "Концентрирование ДНК; запрос сырья",
            )
        elif plant_det == "-" and vpk_det == "+":
            return "ДНК растения не обнаружена (реакция валидна)", ""
        elif plant_det == "мало" and (prom_det == "+" or nos_det == "+"):
            return (
                "ГМ растение; для точного количественного анализа требуется концентрирование ДНК",
                "Концентрирование и идентификация",
            )
        elif plant_det == "-" and vpk_det == "-":
            return "Ложноотрицательный результат", "Перевыделить ДНК и повторить анализ"
        else:
            return "Результат не определён", "Проверить данные"


def analyze_chicken_turkey(protocol: dict, cfg: dict) -> dict:
    samples = protocol["samples"]
    pko_s, oko_s, kov_s = None, None, None
    test_samples = []

    for s in samples:
        n = s["name"].lower().strip()
        if "пко" in n or n.startswith("пко"):
            pko_s = s
        elif "око" in n or n.startswith("око"):
            oko_s = s
        elif "ко-в" in n or "меланж" in n or "10%" in n or "kov" in n:
            kov_s = s
        else:
            test_samples.append(s)

    if not pko_s or not oko_s:
        return {"error": "Не найден ПКО или ОКО. Проверьте названия лунок."}

    CUT = cfg['ct_cutoff']   

    def _det_simple(ct_val) -> str:
        if ct_val is not None and ct_val < CUT:
            return "+"
        return "-"

    def _turkey_text(ct_val) -> str:
        if ct_val is not None and ct_val < CUT:
            return "В пробе содержится ДНК индейки"
        return "В пробе отсутствует ДНК индейки"

    kov_chicken = _ct(kov_s, "Rox") if kov_s else None
    kov_vpk     = _ct(kov_s, "Hex") if kov_s else None
    kov_valid   = (kov_chicken is not None and kov_chicken <= 29) 

    def _chicken_result(s_chicken, s_vpk) -> str:
        if not kov_valid:
            return "Ложноотрицательный результат"
        if s_chicken is not None:
            if s_chicken < kov_chicken:
                return "В пробе содержится ДНК курицы. Количество ДНК соответствует наличию мяса курицы."
            elif s_chicken < CUT:
                return "В пробе содержится ДНК курицы. Количество ДНК менее ДНК 10% меланжа. Мясо курицы (продукт убоя) отсутствует."
            else:
                return "В пробе отсутствует ДНК курицы."
        else:
            return "-"

    pko_turkey  = _ct(pko_s, "Fam")
    pko_chicken = _ct(pko_s, "Rox")
    pko_vpk     = _ct(pko_s, "Hex")

    oko_turkey  = _ct(oko_s, "Fam")
    oko_chicken = _ct(oko_s, "Rox")
    oko_vpk     = _ct(oko_s, "Hex")

    sample_results = []
    for s in test_samples:
        s_turkey  = _ct(s, "Fam")
        s_chicken = _ct(s, "Rox")
        s_vpk     = _ct(s, "Hex")

        turkey_det    = _det_simple(s_turkey)
        turkey_text   = _turkey_text(s_turkey)
        chicken_text  = _chicken_result(s_chicken, s_vpk)
        vpk_det       = _det_simple(s_vpk)

        sample_results.append({
            "name": s["name"], "well": s["well"],
            "turkey_ct": s_turkey,
            "turkey_det": turkey_det,
            "turkey_result": turkey_text,
            "chicken_ct": s_chicken,
            "chicken_result": chicken_text,
            "vpk_ct": s_vpk,
            "vpk_det": vpk_det,
        })

    return {
        "pko": {
            "name": pko_s["name"], "well": pko_s["well"],
            "turkey_ct": pko_turkey,
            "turkey_det": _det_simple(pko_turkey),
            "chicken_ct": pko_chicken,
            "chicken_det": _det_simple(pko_chicken),
            "vpk_ct": pko_vpk,
            "vpk_det": _det_simple(pko_vpk),
        },
        "oko": {
            "name": oko_s["name"], "well": oko_s["well"],
            "turkey_ct": oko_turkey,
            "turkey_det": _det_simple(oko_turkey),
            "chicken_ct": oko_chicken,
            "chicken_det": _det_simple(oko_chicken),
            "vpk_ct": oko_vpk,
            "vpk_det": _det_simple(oko_vpk),
        },
        "kov": {
            "name": kov_s["name"] if kov_s else "КО-В",
            "well": kov_s["well"] if kov_s else "—",
            "chicken_ct": kov_chicken,
            "chicken_det": "+" if kov_valid else "-",
            "vpk_ct": kov_vpk,
            "vpk_det": _det_simple(kov_vpk),
            "valid": kov_valid,
        } if kov_s else None,
        "samples": sample_results,
        "kov_valid": kov_valid,
        "kov_chicken_ct": kov_chicken,
    }


# ═══════════════════════════════════════════════════════════════════
# БЛОК 3: ИНТЕРФЕЙС STREAMLIT
# ═══════════════════════════════════════════════════════════════════

st.markdown("""
<div class="header-card">
    <h1>🧬 RT-PCR Анализатор</h1>
    <p>Автоматическая интерпретация результатов из протоколов амплификатора</p>
</div>
""", unsafe_allow_html=True)

if rtf_to_text is None:
    st.sidebar.warning("⚠️ Библиотека 'striprtf' не найдена! Установите: `pip install striprtf`")

st.sidebar.header("⚙️ Настройки порогов")
test_type_manual = st.sidebar.selectbox(
    "Тест-система:",
    ["Автоопределение", "Растение/35S+FMV/NOS", "Курица/Индейка"]
)

st.sidebar.subheader("🌱 ГМО Скрининг")
is_swab_mode = st.sidebar.checkbox(
    "🧪 Режим смывов (контаминация)",
    help="Отсутствие ДНК растения = норма (смыв чистый)."
)
cfg_gmo = {
    'delta_pos':     st.sidebar.number_input("Δ Ct порог «+» (ПКО +)", value=3.4, step=0.1,
                                              help="Excel: +3.4. Если Ct_образца > Ct_ПКО + X → не детектируется"),
    'delta_low':     st.sidebar.number_input("Δ Ct порог «мало» (ПКО +)", value=10.5, step=0.1,
                                              help="Excel: +10.5. Если Ct > ПКО + X → «-», иначе «мало»"),
    'pko_plant_max': st.sidebar.number_input("Max Ct ПКО (Растение)", value=23.0, step=0.5,
                                              help="Excel: 23. Если Ct_ПКО_plant > X → ПКО не принимается"),
    'pko_prom_max':  st.sidebar.number_input("Max Ct ПКО (Промоторы)", value=36.0, step=0.5,
                                              help="Excel: 36. Если Ct_ПКО_35S/NOS > X → ПКО не принимается"),
    'vpk_max':       st.sidebar.number_input("Max Ct ВПК (ОКО)", value=36.0, step=0.5,
                                              help="Порог ВПК для ОКО"),
    'vpk_delta':     st.sidebar.number_input("Δ Ct ВПК образца к ОКО (+Ct)", value=5.0, step=0.5,
                                              help="Excel: +5. Ct_VPK_образца > Ct_VPK_ОКО + X → «-»"),
}

st.sidebar.subheader("🍗 Курица/Индейка")
cfg_meat = {
    'ct_cutoff': st.sidebar.number_input("Ct cut-off (общий порог)", value=35.0, step=0.5,
                                          help="Excel: <35 → «+», ≥35 → «-»"),
}

# ─── Вкладки ───
tab1, tab2 = st.tabs(["📑 Интерпретация результатов", "📈 Сырые данные и графики"])

protocol_data = None   # глобально передаём во вкладку 2

# ─────────────────────────────────────────
# ВКЛАДКА 1: Интерпретация
# ─────────────────────────────────────────
with tab1:
    st.markdown("### 📑 Загрузка протокола амплификатора")
    uploaded_rtf = st.file_uploader(
        "📂 RTF-протокол (.rtf, .txt)",
        type=["rtf", "txt"],
        key="rtf_uploader"
    )

    if uploaded_rtf:
        try:
            protocol_data = parse_rtf_protocol(uploaded_rtf.read())
        except Exception as e:
            st.error(f"Ошибка чтения файла: {e}")
            protocol_data = None

        if protocol_data is not None:
            if not protocol_data["samples"]:
                st.warning(
                    "В файле не найдены строки с результатами. "
                    "Убедитесь, что таблица содержит заголовок 'Номер лунки' и колонки каналов (Fam, Hex, Rox...)."
                )
            else:
                # ─── Определение типа теста ───
                detected_name = protocol_data["test_name"].lower()
                if test_type_manual != "Автоопределение":
                    active_test = "gmo" if "Растение" in test_type_manual else "meat"
                else:
                    if any(x in detected_name for x in ["растение", "35s", "fmv", "gmo", "nos"]):
                        active_test = "gmo"
                    elif any(x in detected_name for x in ["курица", "индейка", "kurit", "meat"]):
                        active_test = "meat"
                    else:
                        active_test = "unknown"

                # ─── Метаданные протокола ───
                with st.expander("📋 Сведения о протоколе", expanded=True):
                    c1, c2, c3, c4 = st.columns(4)
                    c1.metric("Дата", protocol_data["date"])
                    c2.metric("Оператор", protocol_data["operator"])
                    c3.metric("Протокол №", protocol_data["protocol_number"])
                    c4.metric("Образцов", len(protocol_data["samples"]))
                    st.caption(f"**Тест:** {protocol_data['test_name']}  |  **Программа:** {protocol_data['amplification_program']}")

                st.divider()

                if active_test == "unknown":
                    st.error(
                        "Не удалось автоматически определить тип теста. "
                        "Выберите тест-систему вручную в панели слева."
                    )

                # ══════════════════════════════════════════════
                # РЕЗУЛЬТАТЫ ГМО
                # ══════════════════════════════════════════════
                elif active_test == "gmo":
                    st.markdown("### 🌱 Растение / 35S+FMV / NOS скрининг")
                    if is_swab_mode:
                        st.info("Режим смывов: обнаружение любого сигнала трактуется как контаминация поверхностей.")

                    res = analyze_plant_gmo(protocol_data, cfg_gmo, is_swab=is_swab_mode)

                    if "error" in res:
                        st.error(res["error"])
                    else:
                        # ─── Контроли ───
                        st.markdown("#### Контроли")
                        col1, col2 = st.columns(2)

                        with col1:
                            pko = res["pko"]
                            box = "success-box" if pko["status"] == "Принимается" else "error-box"
                            st.markdown(f"""
                            <div class="{box}">
                                <b>ПКО — {pko['name']}</b> (лунка {pko['well']})<br>
                                Растение (Hex): <code>{pko['plant_ct'] if pko['plant_ct'] else '—'}</code> &nbsp;
                                35S+FMV (Rox): <code>{pko['promoter_ct'] if pko['promoter_ct'] else '—'}</code> &nbsp;
                                NOS (Fam): <code>{pko['nos_ct'] if pko['nos_ct'] else '—'}</code> &nbsp;
                                ВПК (Cy5): <code>{pko['vpk_ct'] if pko['vpk_ct'] else '—'}</code><br>
                                <b>Статус: {pko['status']}</b>{(' — ' + pko['action']) if pko['action'] else ''}
                            </div>
                            """, unsafe_allow_html=True)

                        with col2:
                            oko = res["oko"]
                            box = "success-box" if oko["status"] == "Принимается" else "error-box"
                            st.markdown(f"""
                            <div class="{box}">
                                <b>ОКО — {oko['name']}</b> (лунка {oko['well']})<br>
                                Растение (Hex): <code>{oko['plant_ct'] if oko['plant_ct'] else '—'}</code> &nbsp;
                                ВПК (Cy5): <code>{oko['vpk_ct'] if oko['vpk_ct'] else '—'}</code><br>
                                <b>Статус: {oko['status']}</b>{(' — ' + oko['action']) if oko['action'] else ''}
                            </div>
                            """, unsafe_allow_html=True)

                        st.divider()
                        st.markdown("#### Образцы (Сводная таблица)")

                        rows_gmo = []
                        for sr in res["samples"]:

                            def _fmt(ct, det):
                                if ct:
                                    return f"{ct:.1f} ({det})"
                                return f"— ({det})"

                            rows_gmo.append({
                                "Лунка": sr["well"],
                                "Образец": sr["name"],
                                "Растение (Hex)": _fmt(sr["plant_ct"], sr["plant_det"]),
                                "35S+FMV (Rox)": _fmt(sr["promoter_ct"], sr["promoter_det"]),
                                "NOS (Fam)": _fmt(sr["nos_ct"], sr["nos_det"]),
                                "ВПК (Cy5)": _fmt(sr["vpk_ct"], sr["vpk_det"]),
                                "Вывод": sr["conclusion"],
                                "Действие": sr["action"],
                            })

                        df_gmo = pd.DataFrame(rows_gmo)

                        # Цветовая разметка
                        def _color_gmo(val):
                            if isinstance(val, str):
                                if "(+)" in val:
                                    return "color: #d32f2f; font-weight: bold"
                                if "(-)" in val or "— (-)" in val:
                                    return "color: #388e3c"
                                if "(мало)" in val:
                                    return "color: #f57c00"
                                if "ГМ растение" in val or "КОНТАМИН" in val:
                                    return "background-color: #ffebee; color: #d32f2f; font-weight: bold"
                                if "чистый" in val or "не обнаружена" in val or "не обнаружены" in val:
                                    return "color: #388e3c"
                            return ""

                        if not df_gmo.empty:
                            st.dataframe(
                                df_gmo.style.map(_color_gmo, subset=["Растение (Hex)", "35S+FMV (Rox)", "NOS (Fam)", "ВПК (Cy5)", "Вывод"]),
                                use_container_width=True,
                                hide_index=True,
                            )

                            # Экспорт
                            excel_bytes = to_excel(df_gmo)
                            st.download_button(
                                "⬇️ Скачать результаты (Excel)",
                                data=excel_bytes,
                                file_name=f"GMO_results_{protocol_data['date']}.xlsx",
                                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                            )

                            with st.expander("ℹ️ Биологическая справка: что означают каналы флуоресценции?"):
                                st.markdown("""
                                * 🟡 **HEX (Растение)**: Обнаруживает универсальный ген растения (например, хлоропластную ДНК). Наличие сигнала доказывает, что растительная ДНК успешно выделена из образца и пригодна для ПЦР.
                                * 🟠 **ROX (35S+FMV)**: Выявляет регуляторные элементы — промотор 35S (вирус мозаики цветной капусты) и промотор FMV (вирус мозаики норичника). Эти мощные вирусные "стартеры" искусственно встраивают в геном для запуска работы целевых ГМО-генов.
                                * 🟢 **FAM (NOS)**: Выявляет терминатор NOS (из генома почвенной бактерии *Agrobacterium tumefaciens*). В ГМ-растениях он используется как "стоп-сигнал", завершающий считывание встроенного гена.
                                * 🔴 **Cy5 (ВПК)**: Внутренний положительный контроль. Это искусственная ДНК-мишень, которую добавляют в реакционную смесь. Если ВПК "светится", значит ферменты ПЦР активны и в образце нет ингибиторов (примесей, мешающих реакции).
                                """)

                            st.markdown("#### Подробные выводы")
                            for sr in res["samples"]:
                                if "ОБНАРУЖЕНА КОНТАМИНАЦИЯ" in sr["conclusion"] or "ГМ растение" in sr["conclusion"] or "невалидна" in sr["conclusion"]: 
                                    icon, box = "🔴", "error-box"
                                elif "Смыв чистый" in sr["conclusion"] or "не обнаружена" in sr["conclusion"]: 
                                    icon, box = "🟢", "success-box"
                                elif "Следовая контаминация" in sr["conclusion"] or "Нет растительной" in sr["conclusion"]: 
                                    icon, box = "🟡", "warn-box"
                                else: 
                                    icon, box = "🔵", "info-box"
                                
                                details = []
                                if sr['plant_det'] == '+':
                                    details.append("🟡 <b>HEX (Растение):</b> Обнаружена ДНК растения.")
                                elif sr['plant_det'] == 'мало':
                                    details.append("🟡 <b>HEX (Растение):</b> Обнаружены следы ДНК растения.")
                                
                                if sr['promoter_det'] == '+':
                                    details.append("🟠 <b>ROX (35S+FMV):</b> Выявлены вирусные промоторы (маркер ГМО).")
                                elif sr['promoter_det'] == 'мало':
                                    details.append("🟠 <b>ROX (35S+FMV):</b> Выявлены следы вирусных промоторов (маркер ГМО).")
                                
                                if sr['nos_det'] == '+':
                                    details.append("🟢 <b>FAM (NOS):</b> Выявлен терминатор (маркер ГМО).")
                                elif sr['nos_det'] == 'мало':
                                    details.append("🟢 <b>FAM (NOS):</b> Выявлены следы терминатора (маркер ГМО).")
                                
                                if sr['vpk_det'] == '+':
                                    details.append("🔴 <b>Cy5 (ВПК):</b> ПЦР прошла успешно, ингибиторов нет.")
                                else:
                                    details.append("🔴 <b>Cy5 (ВПК):</b> Сигнал отсутствует (ошибка реакции или ингибирование).")
                                
                                details_html = "<br>".join(details)
                                
                                st.markdown(f"""
                                <div class="{box}">
                                    {icon} <b>{sr['name']}</b> (лунка {sr['well']})<br>
                                    {sr['conclusion']}<br>
                                    <i>Действие: {sr['action']}</i>
                                    <div style="margin-top: 10px; padding-top: 8px; border-top: 1px solid rgba(0,0,0,0.1); font-size: 0.9em; line-height: 1.4;">
                                        <b>Зафиксированные сигналы:</b><br>
                                        {details_html}
                                    </div>
                                </div>
                                """, unsafe_allow_html=True)

                        else:
                            st.info("Нет клинических образцов для отображения (в протоколе найдены только контроли).")

                # ══════════════════════════════════════════════
                # РЕЗУЛЬТАТЫ КУРИЦА/ИНДЕЙКА
                # ══════════════════════════════════════════════
                elif active_test == "meat":
                    st.markdown("### 🍗 Курица / Индейка")
                    res = analyze_chicken_turkey(protocol_data, cfg_meat)

                    if "error" in res:
                        st.error(res["error"])
                    else:
                        # ─── Контроли ───
                        st.markdown("#### Контроли")
                        cols = st.columns(3) if res["kov"] else st.columns(2)

                        with cols[0]:
                            pko = res["pko"]
                            all_ok = all(v == "+" for v in [pko["turkey_det"], pko["chicken_det"], pko["vpk_det"]])
                            box = "success-box" if all_ok else "error-box"
                            st.markdown(f"""
                            <div class="{box}">
                                <b>ПКО — {pko['name']}</b> (лунка {pko['well']})<br>
                                Индейка (Fam): <code>{pko['turkey_ct'] if pko['turkey_ct'] else '—'}</code> → <b>{pko['turkey_det']}</b><br>
                                Курица (Rox): <code>{pko['chicken_ct'] if pko['chicken_ct'] else '—'}</code> → <b>{pko['chicken_det']}</b><br>
                                ВПК (Hex): <code>{pko['vpk_ct'] if pko['vpk_ct'] else '—'}</code> → <b>{pko['vpk_det']}</b>
                            </div>
                            """, unsafe_allow_html=True)

                        with cols[1]:
                            oko = res["oko"]
                            box = "success-box" if oko["turkey_det"] == "-" and oko["chicken_det"] == "-" else "error-box"
                            st.markdown(f"""
                            <div class="{box}">
                                <b>ОКО — {oko['name']}</b> (лунка {oko['well']})<br>
                                Индейка (Fam): <code>{oko['turkey_ct'] if oko['turkey_ct'] else '—'}</code> → <b>{oko['turkey_det']}</b><br>
                                Курица (Rox): <code>{oko['chicken_ct'] if oko['chicken_ct'] else '—'}</code> → <b>{oko['chicken_det']}</b><br>
                                ВПК (Hex): <code>{oko['vpk_ct'] if oko['vpk_ct'] else '—'}</code> → <b>{oko['vpk_det']}</b>
                            </div>
                            """, unsafe_allow_html=True)

                        if res["kov"] and len(cols) > 2:
                            with cols[2]:
                                kov = res["kov"]
                                box = "success-box" if kov["valid"] else "warn-box"
                                st.markdown(f"""
                                <div class="{box}">
                                    <b>КО-В — {kov['name']}</b> (лунка {kov['well']})<br>
                                    Курица (Rox): <code>{kov['chicken_ct'] if kov['chicken_ct'] else '—'}</code> → <b>{kov['chicken_det']}</b><br>
                                    ВПК (Hex): <code>{kov['vpk_ct'] if kov['vpk_ct'] else '—'}</code> → <b>{kov['vpk_det']}</b><br>
                                    {'<b style="color:#388e3c">КО-В принят ✓</b>' if kov['valid'] else '<b style="color:#f57c00">КО-В не принят (Ct &gt; 29)</b>'}
                                </div>
                                """, unsafe_allow_html=True)

                        if not res["kov_valid"]:
                            st.warning(
                                "⚠️ КО-В (меланж 10%) не принят — Ct ROX > 29. "
                                "Все результаты по курице будут отображены как «Ложноотрицательный результат»."
                            )

                        st.divider()
                        st.markdown("#### Образцы (Сводная таблица)")

                        rows_meat = []
                        for sr in res["samples"]:
                            rows_meat.append({
                                "Лунка": sr["well"],
                                "Образец": sr["name"],
                                "Ct Индейка (Fam)": f"{sr['turkey_ct']:.1f}" if sr['turkey_ct'] else "—",
                                "ДНК Индейки": sr["turkey_det"],
                                "Ct Курица (Rox)": f"{sr['chicken_ct']:.1f}" if sr['chicken_ct'] else "—",
                                "Ct ВПК (Hex)": f"{sr['vpk_ct']:.1f}" if sr['vpk_ct'] else "—",
                                "ВПК": sr["vpk_det"],
                                "Заключение по индейке": sr["turkey_result"],
                                "Заключение по курице": sr["chicken_result"],
                            })

                        df_meat = pd.DataFrame(rows_meat)

                        def _color_meat(val):
                            if val == "+":
                                return "color: #d32f2f; font-weight: bold"
                            if val == "-":
                                return "color: #388e3c"
                            if "содержится" in str(val) and "менее" not in str(val):
                                return "color: #d32f2f; font-weight: bold"
                            if "отсутствует" in str(val):
                                return "color: #388e3c"
                            if "Ложноотрицательный" in str(val):
                                return "color: #f57c00; font-weight: bold"
                            return ""

                        if not df_meat.empty:
                            st.dataframe(
                                df_meat.style.map(_color_meat,
                                                  subset=["ДНК Индейки", "ВПК",
                                                          "Заключение по индейке", "Заключение по курице"]),
                                use_container_width=True,
                                hide_index=True,
                            )

                            excel_bytes_m = to_excel(df_meat)
                            st.download_button(
                                "⬇️ Скачать результаты (Excel)",
                                data=excel_bytes_m,
                                file_name=f"Meat_results_{protocol_data['date']}.xlsx",
                                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                            )

                            with st.expander("ℹ️ Биологическая справка: что означают каналы флуоресценции?"):
                                st.markdown("""
                                * 🟢 **FAM (Индейка)**: Обнаруживает строго видоспецифичный генетический маркер, характерный только для индейки (*Meleagris gallopavo*). 
                                * 🟠 **ROX (Курица)**: Обнаруживает строго видоспецифичный генетический маркер домашней курицы (*Gallus gallus*).
                                * 🟡 **HEX (ВПК)**: Внутренний положительный контроль. Это искусственный фрагмент ДНК, необходимый для валидации теста. Мясные и пищевые продукты (особенно со специями) часто содержат вещества, подавляющие ПЦР. Сигнал ВПК подтверждает, что реакция прошла успешно и отсутствие мишеней курицы/индейки — это истинный результат, а не сбой реакции из-за ингибиторов.
                                """)

                            st.markdown("#### Подробные выводы")
                            for sr in res["samples"]:
                                turkey_found = "содержится ДНК индейки" in sr["turkey_result"]
                                chicken_found = "содержится ДНК курицы" in sr["chicken_result"]
                                
                                if turkey_found or chicken_found:
                                    box = "warn-box"
                                    icon = "🔴" if "мяса курицы" in sr["chicken_result"] else "🟡"
                                elif "Ложноотрицательный" in sr["chicken_result"] or "Ошибка" in sr["chicken_result"]:
                                    box = "error-box"
                                    icon = "⚠️"
                                else:
                                    box = "success-box"
                                    icon = "🟢"
                                
                                details = []
                                if turkey_found:
                                    details.append("🟢 <b>FAM (Индейка):</b> Присутствует специфичная ДНК индейки.")
                                if chicken_found:
                                    details.append("🟠 <b>ROX (Курица):</b> Присутствует специфичная ДНК курицы.")
                                
                                if sr['vpk_det'] == '+':
                                    details.append("🟡 <b>HEX (ВПК):</b> ПЦР прошла успешно, ингибиторов нет.")
                                else:
                                    details.append("🟡 <b>HEX (ВПК):</b> Сигнал ВПК отсутствует (ошибка реакции или ингибирование).")
                                
                                details_html = "<br>".join(details)
                                
                                st.markdown(f"""
                                <div class="{box}">
                                    {icon} <b>{sr['name']}</b> (лунка {sr['well']})<br>
                                    Индейка: {sr['turkey_result']}<br>
                                    Курица: {sr['chicken_result']}
                                    <div style="margin-top: 10px; padding-top: 8px; border-top: 1px solid rgba(0,0,0,0.1); font-size: 0.9em; line-height: 1.4;">
                                        <b>Зафиксированные сигналы:</b><br>
                                        {details_html}
                                    </div>
                                </div>
                                """, unsafe_allow_html=True)
                        else:
                            st.info("Нет клинических образцов для отображения (в протоколе найдены только контроли).")

    else:
        # ВОССТАНОВЛЕННЫЕ ИНСТРУКЦИИ НА ПЕРВОЙ ВКЛАДКЕ (При пустом экране)
        st.markdown("---")
        col1, col2 = st.columns(2)
        with col1:
            st.markdown("""
            #### Поддерживаемые тесты
            - **Растение/35S+FMV/NOS скрининг** — скрининг ГМО в растительном сырье
            - **Курица/Индейка** — определение ДНК курицы и индейки в продуктах
            """)
        with col2:
            st.markdown("""
            #### Как пользоваться
            1. Экспортируйте протокол из амплификатора в формате RTF
            2. Загрузите файл через кнопку выше
            3. При необходимости скорректируйте пороги в панели слева
            4. Результаты, контрольные точки и выводы будут рассчитаны автоматически
            """)

# ─────────────────────────────────────────
# ВКЛАДКА 2: Сырые данные и графики
# ─────────────────────────────────────────
with tab2:
    st.markdown("### 📈 Сырые данные амплификации (.r48 / .r96)")
    uploaded_r48 = st.file_uploader(
        "📂 Файл данных (.r48, .r96)",
        type=["r48", "r96"],
        key="r48_uploader"
    )

    if uploaded_r48:
        try:
            raw_dfs = parse_r48(uploaded_r48.read())
        except Exception as e:
            st.error(f"Ошибка чтения файла данных: {e}")
            raw_dfs = {}

        if raw_dfs:
            channels_present = list(raw_dfs.keys())
            st.success(f"Загружено каналов: {', '.join(channels_present)}. "
                       f"Лунок: {len(raw_dfs[channels_present[0]].columns)}.")

            # ─── Сводная таблица Ct ───
            summary_df = calculate_raw_summary(raw_dfs, protocol_data)
            if not summary_df.empty:
                st.markdown("#### Сводная таблица Ct по лункам")
                st.dataframe(summary_df, use_container_width=True, hide_index=True)
                st.download_button(
                    "⬇️ Скачать сводную таблицу (Excel)",
                    data=to_excel(summary_df),
                    file_name="Summary_Ct.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                )

            # ─── Графики кривых амплификации ───
            st.markdown("#### Кривые амплификации")
            
            col1, col2 = st.columns([1, 2])
            with col1:
                view_mode = st.radio("Режим просмотра графиков:", 
                                     ["Один образец (все каналы)", "Один канал (несколько лунок)"],
                                     help="Чтобы увидеть классический мультиплексный график, выберите 'Один образец' и укажите лунку с ПКО.")
                apply_baseline = st.checkbox("Вычесть базовую линию", value=True, help="Строгое вычитание фона по циклам 5-15.")
                apply_smoothing = st.checkbox("Сгладить кривые (Сигмоида)", value=True, help="Применяет математическое сглаживание для устранения 'зубчатости' графика.")
            
            with col2:
                if view_mode == "Один образец (все каналы)":
                    all_wells = raw_dfs[channels_present[0]].columns.tolist()
                    selected_well = st.selectbox("Выберите лунку (рекомендуется ПКО):", all_wells)
                    
                    plot_data = {}
                    for ch in channels_present:
                        if selected_well in raw_dfs[ch].columns:
                            plot_data[ch] = raw_dfs[ch][selected_well]
                    df_display = pd.DataFrame(plot_data)
                else:
                    selected_ch = st.selectbox("Выберите канал флуоресценции:", channels_present)
                    all_wells = raw_dfs[selected_ch].columns.tolist()
                    selected_wells = st.multiselect("Выберите лунки для отображения:", all_wells, default=all_wells[:4])
                    if selected_wells:
                        df_display = raw_dfs[selected_ch][selected_wells]
                    else:
                        df_display = pd.DataFrame()
                    
            if not df_display.empty:
                if apply_smoothing:
                    df_display = df_display.rolling(window=3, min_periods=1, center=True).mean()
                
                if apply_baseline:
                    if len(df_display) > 15:
                        baseline = df_display.loc[5:15].mean()
                    else:
                        baseline = df_display.min()
                    
                    # ПРИНУДИТЕЛЬНОЕ И БЕЗОПАСНОЕ ВЫЧИТАНИЕ (Решает проблему "съехавшего" ROX)
                    df_display = df_display.sub(baseline, axis='columns')

                fig = go.Figure()
                
                # Возвращаем привычные цвета как в бэкапе
                colors_map = {
                    "Fam": "#2ca02c",  # Зеленый
                    "Hex": "#ffb000",  # Желтый
                    "Rox": "#ff5722",  # Оранжевый / Красный
                    "Cy5": "#d32f2f",  # Темно-красный 
                    "Cy5.5": "#7b1fa2" # Фиолетовый
                }

                if view_mode == "Один образец (все каналы)":
                    for ch in df_display.columns:
                        name = selected_well
                        if protocol_data and "samples" in protocol_data:
                            for smp in protocol_data["samples"]:
                                if smp["well"] == selected_well and smp["name"]:
                                    name = f"{selected_well} — {smp['name']}"
                                    break
                        fig.add_trace(go.Scatter(
                            x=df_display.index,
                            y=df_display[ch],
                            name=f"{ch} | {name}",
                            mode='lines',
                            line=dict(color=colors_map.get(ch, "#607D8B"), width=2),
                            opacity=0.85,
                        ))
                else:
                    for w in df_display.columns:
                        name = w
                        if protocol_data and "samples" in protocol_data:
                            for smp in protocol_data["samples"]:
                                if smp["well"] == w and smp["name"]:
                                    name = f"{w} — {smp['name']}"
                                    break
                        fig.add_trace(go.Scatter(
                            x=df_display.index,
                            y=df_display[w],
                            name=f"{selected_ch} | {name}",
                            mode='lines',
                            line=dict(width=2),
                            opacity=0.85,
                        ))

                fig.update_layout(
                    title="Кривые амплификации",
                    xaxis_title="Цикл амплификации",
                    yaxis_title="Флуоресценция (dF)",
                    legend=dict(orientation="v", x=1.02, y=1),
                    height=600,
                    plot_bgcolor='white',
                    hovermode="x unified",
                    xaxis=dict(showgrid=True, gridcolor='#e0e0e0', zeroline=False),
                    yaxis=dict(showgrid=True, gridcolor='#e0e0e0', zeroline=False),
                    margin=dict(l=20, r=180, t=20, b=20)
                )
                st.plotly_chart(fig, use_container_width=True)
        else:
            st.warning("Файл прочитан, но данные каналов не найдены.")
    else:
        if protocol_data:
            st.info(
                "Загрузите файл .r48/.r96 для просмотра кривых амплификации. "
                "RTF-протокол уже загружен — имена образцов будут синхронизированы автоматически."
            )
        else:
            st.info("Загрузите файл .r48 или .r96 для просмотра и анализа сырых кривых амплификации.")