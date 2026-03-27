import streamlit as st
import re
import pandas as pd
import io
import plotly.express as px
import math
from io import BytesIO

try:
    from striprtf.striprtf import rtf_to_text
except ImportError:
    rtf_to_text = None

# Настройка страницы
st.set_page_config(
    page_title="RT-PCR Анализатор",
    page_icon="🧬",
    layout="wide",
)

# ─── Custom CSS (Дизайн и оформление элементов) ───
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

def parse_rtf_protocol(content_bytes: bytes) -> dict:
    """Извлекает данные из текстового RTF-отчета амплификатора."""
    text = None
    for enc in ("cp1251", "utf-8", "latin-1"):
        try:
            raw = content_bytes.decode(enc)
            break
        except UnicodeDecodeError:
            raw = None

    if raw is None:
        raise ValueError("Не удалось декодировать файл. Проверьте кодировку.")

    if rtf_to_text is not None and raw.strip().startswith("{\\rtf"):
        text = rtf_to_text(raw)
    else:
        text = raw

    result = {
        "date": "Неизвестно", "operator": "Неизвестно", "protocol_number": "Неизвестно",
        "result_file": "Неизвестно", "test_name": "Неизвестно", "amplification_program": "Неизвестно",
        "samples": [],
    }

    lines = text.split("\n")
    sample_section = False
    headers_map = {} 
    well_header_idx = -1
    
    for line in lines:
        line_clean = line.strip()
        
        # Парсинг шапки документа (Метаданные теста)
        if "Дата:" in line_clean:
            m = re.search(r'Дата:\|*\|*([^|]+)', line_clean)
            if m: result["date"] = m.group(1).strip()
        elif "Оператор:" in line_clean:
            m = re.search(r'Оператор:\|*\|*([^|]+)', line_clean)
            if m: result["operator"] = m.group(1).strip()
        elif "Номер протокола:" in line_clean:
            m = re.search(r'Номер протокола:\|*\|*([^|]+)', line_clean)
            if m: result["protocol_number"] = m.group(1).strip()
        elif "Тест:" in line_clean:
            m = re.search(r'Тест:\|*\|*([^|]+)', line_clean)
            if m: result["test_name"] = m.group(1).strip()
        elif "Программа амплификации:" in line_clean:
            m = re.search(r'Программа амплификации:\|*\|*([^|]+)', line_clean)
            if m: result["amplification_program"] = m.group(1).strip()

        # Динамический поиск колонок (находит Fam, Hex, Rox независимо от их позиции)
        if "Номер лунки" in line_clean and "Идентификатор" in line_clean:
            sample_section = True
            parts = [p.strip() for p in line_clean.split("|")]
            for idx, p in enumerate(parts):
                p_low = p.lower()
                if "номер лунки" in p_low: well_header_idx = idx
                elif "fam" in p_low or "green" in p_low: headers_map["Fam"] = idx
                elif "hex" in p_low or "yellow" in p_low: headers_map["Hex"] = idx
                elif "rox" in p_low or "orange" in p_low: headers_map["Rox"] = idx
                elif "cy5.5" in p_low: headers_map["Cy5.5"] = idx
                elif "cy5" in p_low or "red" in p_low: headers_map["Cy5"] = idx
            continue

        # Чтение строк с результатами по каждой лунке
        if sample_section and line_clean.startswith("|") and not line_clean.startswith("|*"):
            parts = [p.strip() for p in line_clean.split("|")]
            if len(parts) < 3: continue
            
            well, name = "", ""
            well_idx = -1
            
            for i, p in enumerate(parts):
                if re.match(r'^[A-H]\d{1,2}$', p):
                    well = p
                    well_idx = i
                    break
                    
            if well_idx != -1 and len(parts) > well_idx + 1:
                name_raw = parts[well_idx + 1]
                name = re.sub(r'\s*\(.*?\)\s*$', '', name_raw).strip() 
                
                offset = well_idx - well_header_idx if well_header_idx != -1 else 0
                
                ct_values = {}
                for ch, idx in headers_map.items():
                    actual_idx = idx + offset
                    if actual_idx < len(parts):
                        val = parts[actual_idx]
                        if val and re.match(r'^\d+[,\.]\d+$', val):
                            ct_values[ch] = float(val.replace(",", "."))
                        else:
                            ct_values[ch] = None
                
                if well or name:
                    result["samples"].append({"well": well, "name": name, "ct": ct_values})

        if sample_section and "|*" in line_clean:
            break

    return result

def to_excel(df: pd.DataFrame) -> bytes:
    """Генерация Excel-файла для скачивания результатов."""
    output = BytesIO()
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        df.to_excel(writer, index=False, sheet_name='Результаты')
        worksheet = writer.sheets['Результаты']
        for i, col in enumerate(df.columns):
            column_len = max(df[col].astype(str).map(len).max(), len(col)) + 2
            worksheet.set_column(i, i, column_len)
    return output.getvalue()

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
                if c > 20: num_cycles = c
    if num_cycles == 0:
        num_cycles = 50 
        
    channels = ['Fam', 'Hex', 'Rox', 'Cy5', 'Cy5.5']
    raw_data = {ch: [] for ch in channels}
    
    for line in lines:
        parts = line.strip().split()
        if len(parts) < 10: continue 
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
            # Умный фильтр: отсекаем лишние температурные шаги прибора (берем только элонгацию)
            if len(raw_data[ch]) >= num_cycles * 2:
                reads_per_cycle = len(raw_data[ch]) // num_cycles
                blocks = [raw_data[ch][i::reads_per_cycle][:num_cycles] for i in range(reads_per_cycle)]
                block = max(blocks, key=lambda b: sum(b[-1]))
            else:
                block = raw_data[ch][:num_cycles]
                
            num_wells = len(block[0])
            
            if num_wells == 48: cols, rws = 8, 6
            elif num_wells == 96: cols, rws = 12, 8
            else: cols, rws = num_wells, 1
            
            well_names = []
            for r in range(rws):
                for c in range(cols):
                    well_names.append(f"{chr(65+r)}{c+1}")
            
            if len(well_names) != num_wells:
                well_names = [str(i+1) for i in range(num_wells)]
                
            df = pd.DataFrame(block, columns=well_names)
            df.index = range(1, len(df) + 1)
            df.index.name = 'Цикл'
            result[ch] = df
            
    return result

def calculate_raw_summary(raw_dfs: dict, sample_map: dict) -> pd.DataFrame:
    """
    Создает чистую итоговую таблицу Ct из сырых данных.
    Использует классический метод базовой линии (с 5 по 15 цикл), что
    полностью решает проблему ложных "U-образных" кривых.
    """
    if not raw_dfs: return pd.DataFrame()
    first_ch = list(raw_dfs.keys())[0]
    wells = raw_dfs[first_ch].columns.tolist()
    
    rows = []
    display_channels = ["Fam", "Hex", "Rox", "Cy5", "Cy5.5"]
    
    for w in wells:
        sample_name = sample_map.get(w, "")
        row = {"Лунка": w, "Образец": sample_name}
        has_signal = False
        
        for ch in display_channels:
            ct_val = "—"
            if ch in raw_dfs and w in raw_dfs[ch].columns:
                series = raw_dfs[ch][w]
                
                # 1. Сглаживаем кривую (скользящее среднее)
                smooth = series.rolling(window=3, min_periods=1, center=True).mean()
                
                # 2. Вычисляем классическую базовую линию (среднее значение с 5 по 15 цикл)
                # Это идеально защищает от артефактов (провалов и ям)
                if len(smooth) > 15:
                    baseline = smooth.loc[5:15].mean()
                else:
                    baseline = smooth.min()
                    
                norm = smooth - baseline
                
                max_peak = norm.max()
                peak_idx = norm.idxmax()
                
                # УМНЫЕ ФИЛЬТРЫ S-ОБРАЗНОЙ КРИВОЙ ПЦР:
                # Правило 1: Прирост минимум 10 единиц (поймает слабый ROX в A6, но отсечет фоновый шум)
                # Правило 2: Пик строго во второй половине реакции
                if max_peak >= 10.0 and peak_idx > 15:
                    
                    # Правило 3: Производная (Крутизна)
                    # Настоящая ПЦР растет экспоненциально. 
                    # Линейный дрейф пластика будет забракован (max_slope < 1.5)
                    max_slope = norm.diff(periods=2).max()
                    
                    # Правило 4: В конце график не должен рухнуть ниже 5 единиц
                    final_val = norm.iloc[-3:].mean()
                    
                    if max_slope >= 1.5 and final_val >= 5.0:
                        has_signal = True
                        
                        # Вычисляем Ct (пересечение 10% порога, но не ниже 5 RFU)
                        threshold = max(max_peak * 0.1, 5.0)
                        cross = norm[(norm.index <= peak_idx) & (norm >= threshold)]
                        
                        # Защита: Ct не может стоять раньше 5 цикла
                        if not cross.empty and cross.index[0] > 5:
                            ct_val = round(cross.index[0], 1)
                            
            row[f"{ch} Ct"] = ct_val
            
        # Заносим лунку в таблицу, если она не пустая (имеет образец или любой подтвержденный сигнал)
        if sample_name or has_signal:
            if not row["Образец"]:
                row["Образец"] = "Неизвестный образец"
            rows.append(row)
            
    return pd.DataFrame(rows)


# ═══════════════════════════════════════════════════════════════════
# БЛОК 2: АНАЛИЗАТОРЫ ЛОГИКИ ТЕСТ-СИСТЕМ (МАТЕМАТИКА EXCEL)
# ═══════════════════════════════════════════════════════════════════

def _get_ct(sample: dict, ch_name: str) -> float | None:
    """Безопасное извлечение Ct значения из словаря лунки."""
    if not sample: return None
    return sample.get("ct", {}).get(ch_name)

def analyze_plant_gmo(protocol: dict, cfg: dict, is_swab: bool = False) -> dict:
    """Интерпретация скрининга ГМО (Растение, 35S+FMV, NOS)."""
    samples = protocol["samples"]
    pko, oko = None, None
    test_samples = []
    
    # Поиск контролей в списке лунок
    for s in samples:
        n = s["name"].lower()
        if "пко" in n: pko = s
        elif "око" in n: oko = s
        else: test_samples.append(s)
        
    if not pko or not oko:
        return {"error": "Отсутствует ПКО или ОКО. Анализ невозможен. Проверьте названия в протоколе."}

    # Анализ ПКО (Положительный контроль)
    pko_plant = _get_ct(pko, "Hex")
    pko_35s = _get_ct(pko, "Rox")
    pko_nos = _get_ct(pko, "Fam")
    pko_vpk = _get_ct(pko, "Cy5")
    
    pko_status = "Принимается"
    pko_action = ""
    if any(v is None for v in [pko_plant, pko_35s, pko_nos]) or \
       (pko_plant > cfg['pko_plant_max']) or (pko_35s > cfg['pko_prom_max']) or (pko_nos > cfg['pko_prom_max']):
        pko_status = "Не принимается ПКО"
        pko_action = "Повторить анализ"

    # Анализ ОКО (Отрицательный контроль)
    oko_plant = _get_ct(oko, "Hex")
    oko_35s = _get_ct(oko, "Rox")
    oko_nos = _get_ct(oko, "Fam")
    oko_vpk = _get_ct(oko, "Cy5")
    
    oko_has_signal = False
    if oko_plant is not None and (pko_plant is None or oko_plant <= pko_plant + cfg['delta_low']): oko_has_signal = True
    if oko_35s is not None and (pko_35s is None or oko_35s <= pko_35s + cfg['delta_pos']): oko_has_signal = True
    if oko_nos is not None and (pko_nos is None or oko_nos <= pko_nos + cfg['delta_pos']): oko_has_signal = True
    
    vpk_limit = (pko_vpk + cfg['vpk_delta']) if pko_vpk else cfg['vpk_max']
    oko_vpk_ok = (oko_vpk is not None and oko_vpk <= vpk_limit)
    
    oko_status = "Принимается"
    oko_action = ""
    if oko_has_signal: 
        oko_status = "Контаминация"
        oko_action = "Провести деконтаминацию"
    elif not oko_vpk_ok: 
        oko_status = "Не принимается ВПК"
        oko_action = "Повторить анализ"

    # Анализ клинических образцов
    sample_results = []
    for s in test_samples:
        s_plant = _get_ct(s, "Hex")
        s_35s = _get_ct(s, "Rox")
        s_nos = _get_ct(s, "Fam")
        s_vpk = _get_ct(s, "Cy5")
        
        plant_det, prom_det, nos_det = "-", "-", "-"
        
        if is_swab:
            # ЛОГИКА ДЛЯ СМЫВОВ: Любой подтвержденный сигнал = контаминация
            if s_plant is not None: plant_det = "+" if s_plant <= 30 else "мало"
            if s_35s is not None: prom_det = "+" if s_35s <= 30 else "мало"
            if s_nos is not None: nos_det = "+" if s_nos <= 30 else "мало"
        else:
            # ЛОГИКА ДЛЯ ПРОДУКТОВ ПИТАНИЯ (Сравнение дельт с ПКО)
            if s_plant is not None and pko_plant is not None:
                if s_plant <= pko_plant + cfg['delta_pos']: plant_det = "+"
                elif s_plant <= pko_plant + cfg['delta_low']: plant_det = "мало"
                
            if plant_det == "+":
                if s_35s and pko_35s and pko_plant and s_plant:
                    if (s_35s - s_plant) <= (pko_35s - pko_plant) + cfg['delta_pos']: prom_det = "+"
                if s_nos and pko_nos and pko_plant and s_plant:
                    if (s_nos - s_plant) <= (pko_nos - pko_plant) + cfg['delta_pos']: nos_det = "+"
            elif plant_det == "мало":
                if s_35s and pko_35s and s_35s <= pko_35s + cfg['delta_pos']: prom_det = "+"
                if s_nos and pko_nos and s_nos <= pko_nos + cfg['delta_pos']: nos_det = "+"

        # Валидация ВПК
        vpk_det = "+" if (s_vpk and oko_vpk and s_vpk <= oko_vpk + cfg['vpk_delta']) or (s_vpk and not oko_vpk and s_vpk <= cfg['vpk_max']) else "-"
        
        conclusion, action = "Результат не определён", "Проверить данные"
        
        # Формирование итоговых выводов (Тексты из Excel Лист2)
        if is_swab:
            if plant_det == "+" or prom_det == "+" or nos_det == "+":
                conclusion, action = "ОБНАРУЖЕНА КОНТАМИНАЦИЯ!", "Провести повторную деконтаминацию"
            elif plant_det == "мало" or prom_det == "мало" or nos_det == "мало":
                conclusion, action = "Следовая контаминация (мало ДНК)", "Обратить внимание / повторная деконтаминация"
            elif plant_det == "-" and prom_det == "-" and nos_det == "-" and vpk_det == "+":
                conclusion, action = "Смыв чистый (контаминации нет)", "Не требуется"
            elif vpk_det == "-":
                conclusion, action = "Реакция невалидна (ВПК не прошел)", "Повторить анализ"
        else:
            if plant_det == "+" and (prom_det == "+" or nos_det == "+"):
                conclusion, action = "ГМ растение; количество ДНК выше предела количественного определения", "Идентификация и определение количества %"
            elif plant_det == "+" and prom_det == "-" and nos_det == "-":
                conclusion, action = "Промоторы 35S, FMV и терминатор NOS не обнаружены", "Обнаружение сои линий BPS-CV-127-9 и MON87701"
            elif plant_det == "мало" and prom_det == "-" and nos_det == "-":
                conclusion, action = "Промоторы 35S, FMV и терминатор NOS не обнаружены", "Концентрирование ДНК; запрос сырья"
            elif plant_det == "-" and vpk_det == "+":
                conclusion, action = "ДНК растения не обнаружена (реакция валидна)", ""
            elif plant_det == "мало" and (prom_det == "+" or nos_det == "+"):
                conclusion, action = "ГМ растение; для точного количественного анализа требуется концентрирование ДНК", "Концентрирование и идентификация"
            elif plant_det == "-" and prom_det == "-" and nos_det == "-" and vpk_det == "-":
                conclusion, action = "Нет растительной ДНК (Ложноотрицательный)", "Перевыделить ДНК"

        sample_results.append({
            "name": s["name"], "well": s["well"],
            "plant_ct": s_plant, "plant_det": plant_det,
            "promoter_ct": s_35s, "promoter_det": prom_det,
            "nos_ct": s_nos, "nos_det": nos_det,
            "vpk_ct": s_vpk, "vpk_det": vpk_det,
            "conclusion": conclusion, "action": action
        })

    return {
        "pko": {
            "name": pko["name"], "well": pko["well"],
            "plant_ct": pko_plant, "promoter_ct": pko_35s, "nos_ct": pko_nos, "vpk_ct": pko_vpk,
            "status": pko_status, "action": pko_action
        },
        "oko": {
            "name": oko["name"], "well": oko["well"], "vpk_ct": oko_vpk,
            "status": oko_status, "action": oko_action
        },
        "samples": sample_results
    }


def analyze_chicken_turkey(protocol: dict, cfg: dict) -> dict:
    """Интерпретация качественного анализа Курица/Индейка."""
    samples = protocol["samples"]
    pko, oko, kov = None, None, None
    test_samples = []
    
    for s in samples:
        n = s["name"].lower()
        if "пко" in n: pko = s
        elif "око" in n: oko = s
        elif "ко-в" in n or "меланж" in n or "10%" in n: kov = s
        elif "поверка" in n: pass # Пропускаем технологические тесты амплификатора
        else: test_samples.append(s)

    if not pko or not oko:
        return {"error": "Отсутствует ПКО или ОКО. Проверьте названия лунок."}
        
    pko_turkey  = _get_ct(pko, "Fam")
    pko_chicken = _get_ct(pko, "Rox")
    pko_vpk     = _get_ct(pko, "Hex")
    
    oko_turkey  = _get_ct(oko, "Fam")
    oko_chicken = _get_ct(oko, "Rox")
    oko_vpk     = _get_ct(oko, "Hex")
    
    kov_chicken = _get_ct(kov, "Rox") if kov else None
    kov_vpk     = _get_ct(kov, "Hex") if kov else None

    sample_results = []
    for s in test_samples:
        s_turkey = _get_ct(s, "Fam")
        s_chicken = _get_ct(s, "Rox")
        s_vpk = _get_ct(s, "Hex")
        
        # Точное соответствие логике Excel для Индейки
        if s_turkey is not None and s_turkey < cfg['ct_cutoff']:
            turkey_res = "В пробе содержится ДНК индейки"
        else:
            turkey_res = "В пробе отсутствует ДНК индейки"
            
        # Точное соответствие логике Excel для Курицы
        if kov_chicken is None or kov_chicken > 29:
             chicken_res = "Ошибка калибратора (КО-В)"
        elif s_chicken is not None and s_chicken < cfg['ct_cutoff']:
             if s_chicken < kov_chicken:
                 chicken_res = "В пробе содержится ДНК курицы. Количество ДНК соответствует наличию мяса курицы."
             else:
                 chicken_res = "В пробе содержится ДНК курицы. Количество ДНК менее ДНК 10% меланжа. Мясо курицы (продукт убоя) отсутствует."
        else:
             if s_vpk is not None and s_vpk < cfg['ct_cutoff']:
                 chicken_res = "В пробе отсутствует ДНК курицы."
             else:
                 chicken_res = "Ложноотрицательный результат"

        vpk_det = "+" if (s_vpk is not None and s_vpk < cfg['ct_cutoff']) else "-"

        sample_results.append({
            "name": s["name"], "well": s["well"],
            "turkey_ct": s_turkey, "turkey_result": turkey_res,
            "chicken_ct": s_chicken, "chicken_result": chicken_res,
            "vpk_ct": s_vpk, "vpk_det": vpk_det
        })

    return {
        "pko": {
            "name": pko["name"], "well": pko["well"],
            "turkey_ct": pko_turkey, "chicken_ct": pko_chicken, "vpk_ct": pko_vpk,
            "turkey_det": "+" if pko_turkey and pko_turkey < cfg['ct_cutoff'] else "-",
            "chicken_det": "+" if pko_chicken and pko_chicken < cfg['ct_cutoff'] else "-",
            "vpk_det": "+" if pko_vpk and pko_vpk < cfg['ct_cutoff'] else "-"
        },
        "oko": {
            "name": oko["name"], "well": oko["well"],
            "turkey_ct": oko_turkey, "chicken_ct": oko_chicken, "vpk_ct": oko_vpk,
            "turkey_det": "+" if oko_turkey and oko_turkey < cfg['ct_cutoff'] else "-",
            "chicken_det": "+" if oko_chicken and oko_chicken < cfg['ct_cutoff'] else "-",
            "vpk_det": "+" if oko_vpk and oko_vpk < cfg['ct_cutoff'] else "-"
        },
        "kov": {
            "name": kov["name"] if kov else "КО-В", "well": kov["well"] if kov else "-",
            "chicken_ct": kov_chicken, "vpk_ct": kov_vpk,
            "chicken_det": "-" if kov_chicken and kov_chicken > 29 else "+",
            "vpk_det": "+" if kov_vpk and kov_vpk < cfg['ct_cutoff'] else "-"
        } if kov else None,
        "samples": sample_results
    }


# ═══════════════════════════════════════════════════════════════════
# БЛОК 3: ИНТЕРФЕЙС ПРИЛОЖЕНИЯ STREAMLIT
# ═══════════════════════════════════════════════════════════════════

st.markdown("""
<div class="header-card">
    <h1>🧬 RT-PCR Анализатор</h1>
    <p>Автоматическая интерпретация результатов из протоколов амплификатора</p>
</div>
""", unsafe_allow_html=True)

# Боковая панель для настроек алгоритма
st.sidebar.header("⚙️ Настройки порогов")
test_type_manual = st.sidebar.selectbox("Выберите тест-систему:", ["Автоопределение", "Растение/35S+FMV/NOS", "Курица/Индейка"])

st.sidebar.subheader("🌱 ГМО Скрининг")
is_swab_mode = st.sidebar.checkbox("🧪 Режим смывов (на контаминацию)", help="Отсутствие ДНК растения считается нормой (смыв чистый).")
cfg_gmo = {
    'delta_pos': st.sidebar.number_input("Δ Ct Положительный (ПКО+)", value=3.4, step=0.1, help="Порог дельты для детекции"),
    'delta_low': st.sidebar.number_input("Δ Ct Порог 'Мало'", value=10.5, step=0.1),
    'pko_plant_max': st.sidebar.number_input("Max Ct ПКО (Растение)", value=23.0, step=0.5),
    'pko_prom_max': st.sidebar.number_input("Max Ct ПКО (Промоторы)", value=36.0, step=0.5),
    'vpk_max': st.sidebar.number_input("Max Ct ВПК", value=36.0, step=0.5),
    'vpk_delta': st.sidebar.number_input("Допуск ВПК к ОКО (+Ct)", value=5.0, step=0.5),
}

st.sidebar.subheader("🍗 Курица/Индейка")
cfg_meat = {
    'ct_cutoff': st.sidebar.number_input("Порог отсечения (Ct Cut-off)", value=35.0, step=0.5),
}

# Организация вкладок
tab1, tab2 = st.tabs(["📑 Интерпретация результатов", "📈 Сырые данные и графики"])

# Глобальная переменная для передачи имен образцов из 1-й вкладки во 2-ю
protocol_data = None

# ----------------- ВКЛАДКА 1: RTF ИНТЕРПРЕТАЦИЯ -----------------
with tab1:
    st.markdown("### 📑 Загрузка протокола для интерпретации")
    uploaded = st.file_uploader("📂 Загрузите протокол амплификатора (.rtf, .txt)", type=["rtf", "txt"], key="rtf_uploader")

    if uploaded:
        try:
            protocol_data = parse_rtf_protocol(uploaded.read())
        except Exception as e:
            st.error(f"Ошибка чтения файла: {e}")
            protocol_data = None
            
        if protocol_data and not protocol_data["samples"]:
            st.warning("В файле не найдены результаты образцов. Убедитесь, что таблица содержит колонки Fam, Hex, Rox и т.д.")
        elif protocol_data:
            detected_name = protocol_data["test_name"].lower()
            if test_type_manual != "Автоопределение":
                active_test = "gmo" if "Растение" in test_type_manual else "meat"
            else:
                active_test = "gmo" if any(x in detected_name for x in ["растение", "35s", "fmv", "gmo"]) else "meat" if any(x in detected_name for x in ["курица", "индейка", "kurit"]) else "unknown"

            with st.expander("📋 Сведения о протоколе", expanded=True):
                c1, c2, c3 = st.columns(3)
                c1.write(f"**Дата:** {protocol_data['date']}")
                c2.write(f"**Оператор:** {protocol_data['operator']}")
                c3.write(f"**Тест в файле:** {protocol_data['test_name']}")

            if active_test == "unknown":
                st.error("Не удалось определить тип набора. Выберите его вручную в панели слева.")
            else:
                st.divider()
                df_export = pd.DataFrame()
                
                # ==========================================
                # Отрисовка результатов ГМО
                # ==========================================
                if active_test == "gmo":
                    st.markdown("### 🌱 Результаты скрининга ГМО")
                    res = analyze_plant_gmo(protocol_data, cfg_gmo, is_swab=is_swab_mode)
                    if "error" in res:
                        st.error(res["error"])
                    else:
                        st.markdown("#### Контроли")
                        c1, c2 = st.columns(2)
                        with c1:
                            pko = res["pko"]
                            box_class = "success-box" if "Принимается" in pko["status"] else "error-box"
                            st.markdown(f"""
                            <div class="{box_class}">
                                <b>ПКО ({pko['name']})</b> — лунка {pko['well']}<br>
                                Растение: <code>{pko['plant_ct'] or '—'}</code> | 35S+FMV: <code>{pko['promoter_ct'] or '—'}</code> | NOS: <code>{pko['nos_ct'] or '—'}</code> | ВПК: <code>{pko['vpk_ct'] or '—'}</code><br>
                                <b>Статус: {pko['status']}</b> {(' — ' + pko['action']) if pko['action'] else ''}
                            </div>
                            """, unsafe_allow_html=True)
                        with c2:
                            oko = res["oko"]
                            box_class = "success-box" if "Принимается" in oko["status"] else "error-box"
                            st.markdown(f"""
                            <div class="{box_class}">
                                <b>ОКО ({oko['name']})</b> — лунка {oko['well']}<br>
                                ВПК: <code>{oko['vpk_ct'] or '—'}</code><br>
                                <b>Статус: {oko['status']}</b> {(' — ' + oko['action']) if oko['action'] else ''}
                            </div>
                            """, unsafe_allow_html=True)
                            
                        st.markdown("#### Образцы (Сводная таблица)")
                        rows = []
                        for sr in res["samples"]:
                            rows.append({
                                "Лунка": sr["well"], "Образец": sr["name"],
                                "Растение Ct": f"{sr['plant_ct']:.1f}" if sr['plant_ct'] else "—", "Растение": sr["plant_det"],
                                "35S+FMV Ct": f"{sr['promoter_ct']:.1f}" if sr['promoter_ct'] else "—", "35S+FMV": sr["promoter_det"],
                                "NOS Ct": f"{sr['nos_ct']:.1f}" if sr['nos_ct'] else "—", "NOS": sr["nos_det"],
                                "ВПК Ct": f"{sr['vpk_ct']:.1f}" if sr['vpk_ct'] else "—", "ВПК": sr["vpk_det"],
                                "Вывод": sr["conclusion"], "Действие": sr["action"],
                            })
                        df_export = pd.DataFrame(rows)
                        st.dataframe(df_export.style.applymap(lambda x: "background-color: #ffcdd2; font-weight: bold" if x=="+" else ("background-color: #fff9c4; font-weight: bold" if x=="мало" else ""), subset=['Растение', '35S+FMV', 'NOS', 'ВПК']), use_container_width=True)

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

                # ==========================================
                # Отрисовка результатов Курица/Индейка
                # ==========================================
                elif active_test == "meat":
                    st.markdown("### 🍗 Результаты определения курицы / индейки")
                    res = analyze_chicken_turkey(protocol_data, cfg_meat)
                    if "error" in res:
                        st.error(res["error"])
                    else:
                        st.markdown("#### Контроли")
                        cols = st.columns(3 if res.get("kov") else 2)
                        with cols[0]:
                            pko = res["pko"]
                            st.markdown(f"""
                            <div class="info-box">
                                <b>ПКО ({pko['name']})</b> — {pko['well']}<br>
                                Индейка (FAM): <code>{pko['turkey_ct'] or '—'}</code> [{pko['turkey_det']}]<br>
                                Курица (ROX): <code>{pko['chicken_ct'] or '—'}</code> [{pko['chicken_det']}]<br>
                                ВПК (HEX): <code>{pko['vpk_ct'] or '—'}</code> [{pko['vpk_det']}]
                            </div>
                            """, unsafe_allow_html=True)
                        
                        if res.get("kov"):
                            with cols[1]:
                                kov = res["kov"]
                                st.markdown(f"""
                                <div class="{'success-box' if kov['chicken_det'] == '+' else 'error-box'}">
                                    <b>КО-В ({kov['name']})</b> — {kov['well']}<br>
                                    Курица (ROX): <code>{kov['chicken_ct'] or '—'}</code> [{kov['chicken_det']}]<br>
                                    ВПК (HEX): <code>{kov['vpk_ct'] or '—'}</code> [{kov['vpk_det']}]
                                </div>
                                """, unsafe_allow_html=True)
                        
                        with cols[-1]:
                            oko = res["oko"]
                            st.markdown(f"""
                            <div class="info-box">
                                <b>ОКО ({oko['name']})</b> — {oko['well']}<br>
                                Индейка: <code>{oko['turkey_ct'] or '—'}</code> [{oko['turkey_det']}]<br>
                                Курица: <code>{oko['chicken_ct'] or '—'}</code> [{oko['chicken_det']}]<br>
                                ВПК: <code>{oko['vpk_ct'] or '—'}</code> [{oko['vpk_det']}]
                            </div>
                            """, unsafe_allow_html=True)

                        st.markdown("#### Образцы (Сводная таблица)")
                        rows = []
                        for sr in res["samples"]:
                            rows.append({
                                "Лунка": sr["well"], "Образец": sr["name"],
                                "Индейка Ct (FAM)": f"{sr['turkey_ct']:.1f}" if sr['turkey_ct'] else "—", "Индейка": sr["turkey_result"],
                                "Курица Ct (ROX)": f"{sr['chicken_ct']:.1f}" if sr['chicken_ct'] else "—", "Курица": sr["chicken_result"],
                                "ВПК Ct (HEX)": f"{sr['vpk_ct']:.1f}" if sr['vpk_ct'] else "—", "ВПК": sr["vpk_det"],
                            })
                        df_export = pd.DataFrame(rows)
                        st.dataframe(df_export, use_container_width=True)

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

                if not df_export.empty:
                    st.markdown("<br>", unsafe_allow_html=True)
                    excel_data = to_excel(df_export)
                    st.download_button(
                        label="💾 Скачать отчет (Excel)",
                        data=excel_data,
                        file_name=f"Отчет_{protocol_data['date'].replace(' ', '_')}.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )

    else:
        st.markdown("---")
        col1, col2 = st.columns(2)
        with col1:
            st.markdown("""
            #### Поддерживаемые тесты
            - **Растение/35S+FMV/NOS скрининг** — скрининг ГМО в растительном сырье
            - **Курица/Индейка** — определение ДНК курицы и индейки в продуктов
            """)
        with col2:
            st.markdown("""
            #### Как пользоваться
            1. Экспортируйте протокол из амплификатора в формате RTF
            2. Загрузите файл через кнопку выше
            3. При необходимости скорректируйте пороги в панели слева
            4. Результаты, контрольные точки и выводы будут рассчитаны автоматически
            """)

# ----------------- ВКЛАДКА 2: ГРАФИКИ И СЫРЫЕ ДАННЫЕ -----------------
with tab2:
    st.markdown("### 📈 Визуализация кривых амплификации")
    st.info("Загрузите файл с сырыми данными (например, **.r48** или **.r96** из ПО ДНК-Технология), чтобы посмотреть графики флуоресценции.")
    
    uploaded_raw = st.file_uploader("📂 Загрузите файл сырых данных (.r48, .r96)", type=["r48", "r96", "txt"], key="r48_uploader")
    
    if uploaded_raw:
        try:
            raw_dfs = parse_r48(uploaded_raw.read())
        except Exception as e:
            st.error(f"Ошибка чтения файла данных: {e}")
            raw_dfs = None
            
        if raw_dfs:
            available_channels = list(raw_dfs.keys())
            if not available_channels:
                st.warning("Не удалось найти данные амплификации в файле.")
            else:
                col1, col2 = st.columns([1, 2])
                with col1:
                    view_mode = st.radio("Режим просмотра графиков:", 
                                         ["Один образец (все каналы)", "Один канал (несколько лунок)"],
                                         help="Чтобы увидеть классический мультиплексный график (как в ПО амплификатора), выберите 'Один образец' и укажите лунку с ПКО.")
                    apply_baseline = st.checkbox("Вычесть базовую линию", value=True, help="Строгое вычитание фона по циклам 5-15 (как в ПО амплификатора). Устраняет ямы и дрейф.")
                    apply_smoothing = st.checkbox("Сгладить кривые (Сигмоида)", value=True, help="Применяет математическое сглаживание для устранения 'зубчатости' графика.")
                
                with col2:
                    if view_mode == "Один образец (все каналы)":
                        all_wells = raw_dfs[available_channels[0]].columns.tolist()
                        selected_well = st.selectbox("Выберите лунку (рекомендуется ПКО):", all_wells)
                        
                        plot_data = {}
                        for ch in available_channels:
                            if selected_well in raw_dfs[ch].columns:
                                plot_data[ch] = raw_dfs[ch][selected_well]
                        df_display = pd.DataFrame(plot_data)
                    else:
                        selected_ch = st.selectbox("Выберите канал флуоресценции:", available_channels)
                        all_wells = raw_dfs[selected_ch].columns.tolist()
                        selected_wells = st.multiselect("Выберите лунки для отображения:", all_wells, default=all_wells[:4])
                        df_display = raw_dfs[selected_ch][selected_wells] if selected_wells else pd.DataFrame()
                        
                if not df_display.empty:
                    # Математика графика полностью синхронизирована с математикой таблицы!
                    
                    # 1. Сглаживание
                    if apply_smoothing:
                        df_display = df_display.rolling(window=3, min_periods=1, center=True).mean()
                        
                    # 2. Вычитание базовой линии (циклы с 5 по 15)
                    if apply_baseline:
                        if len(df_display) > 15:
                            baseline = df_display.loc[5:15].mean()
                        else:
                            baseline = df_display.min()
                        df_display = df_display - baseline
                        
                    # 3. Отрисовка
                    if view_mode == "Один образец (все каналы)":
                        color_map = {
                            "Fam": "#2ca02c",  # Зеленый
                            "Hex": "#ffb000",  # Желтый / Оранжево-желтый
                            "Rox": "#ff5722",  # Красный / Глубокий оранжевый
                            "Cy5": "#d32f2f",  # Темно-красный
                            "Cy5.5": "#7b1fa2" # Фиолетовый
                        }
                        fig = px.line(df_display, x=df_display.index, y=df_display.columns, 
                                      color_discrete_map=color_map,
                                      labels={"value": "Уровень флуоресценции (dF)", "Цикл": "Цикл амплификации", "variable": "Канал"})
                    else:
                        fig = px.line(df_display, x=df_display.index, y=df_display.columns,
                                      labels={"value": "Уровень флуоресценции (dF)", "Цикл": "Цикл амплификации", "variable": "Лунка"})
                    
                    fig.update_layout(
                        plot_bgcolor='white',
                        hovermode="x unified",
                        legend_title_text="Легенда",
                        xaxis=dict(showgrid=True, gridcolor='#e0e0e0', zeroline=False),
                        yaxis=dict(showgrid=True, gridcolor='#e0e0e0', zeroline=False),
                        margin=dict(l=20, r=20, t=20, b=20)
                    )
                    
                    st.plotly_chart(fig, use_container_width=True)
                    
                # ── ИТОГОВАЯ ТАБЛИЦА (Из сырых данных) ──
                st.markdown("---")
                st.markdown("#### 📊 Итоговая таблица значений (из сырых данных)")
                
                sample_map = {}
                if protocol_data and "samples" in protocol_data:
                    st.info("💡 Названия образцов автоматически подтянуты из загруженного RTF-протокола на первой вкладке.")
                    for s in protocol_data["samples"]:
                        sample_map[s["well"]] = s["name"]
                else:
                    st.info("💡 Загрузите RTF-протокол на первой вкладке, чтобы в этой таблице отображались названия образцов, а не только номера лунок.")
                        
                summary_df = calculate_raw_summary(raw_dfs, sample_map)
                
                if not summary_df.empty:
                    st.dataframe(summary_df, use_container_width=True)
                    
                    excel_raw = to_excel(summary_df)
                    st.download_button(
                        label="💾 Скачать таблицу сырых данных (Excel)",
                        data=excel_raw,
                        file_name="Сырые_данные_Ct.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )