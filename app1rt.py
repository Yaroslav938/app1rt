import streamlit as st
import re
import pandas as pd
import io

try:
    from striprtf.striprtf import rtf_to_text
except ImportError:
    rtf_to_text = None

st.set_page_config(
    page_title="RT-PCR Анализатор",
    page_icon="🧬",
    layout="wide",
)

# ─── Custom CSS (Возвращены все стили из app1rt.py) ───
st.markdown("""
<style>
    .main .block-container { max-width: 1200px; padding-top: 2rem; }
    .stTabs [data-baseweb="tab-list"] { gap: 8px; }
    .stTabs [data-baseweb="tab"] { padding: 10px 24px; font-weight: 600; }
    .result-positive { color: #d32f2f; font-weight: bold; }
    .result-negative { color: #388e3c; font-weight: bold; }
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
# УТИЛИТЫ (Оптимизированный надежный парсер и Экспорт)
# ═══════════════════════════════════════════════════════════════════

def parse_rtf_protocol(content_bytes: bytes) -> dict:
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
    
    for line in lines:
        line_clean = line.strip()
        
        # Метаданные
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

        # Поиск заголовка таблицы (динамический поиск каналов - это защита от багов экспорта)
        if "Номер лунки" in line_clean and "Идентификатор" in line_clean:
            sample_section = True
            parts = [p.strip() for p in line_clean.split("|")]
            for idx, p in enumerate(parts):
                p_low = p.lower()
                if "fam" in p_low or "green" in p_low: headers_map["Fam"] = idx
                elif "hex" in p_low or "yellow" in p_low: headers_map["Hex"] = idx
                elif "rox" in p_low or "orange" in p_low: headers_map["Rox"] = idx
                elif "cy5.5" in p_low: headers_map["Cy5.5"] = idx
                elif "cy5" in p_low or "red" in p_low: headers_map["Cy5"] = idx
            continue

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
                
                ct_values = {}
                for ch, idx in headers_map.items():
                    if idx < len(parts):
                        val = parts[idx]
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
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        df.to_excel(writer, index=False, sheet_name='Результаты')
        worksheet = writer.sheets['Результаты']
        for i, col in enumerate(df.columns):
            column_len = max(df[col].astype(str).map(len).max(), len(col)) + 2
            worksheet.set_column(i, i, column_len)
    return output.getvalue()


# ═══════════════════════════════════════════════════════════════════
# АНАЛИЗАТОРЫ (С детальным возвратом контролей как в app1rt.py)
# ═══════════════════════════════════════════════════════════════════

def _get_ct(sample: dict, ch_name: str) -> float | None:
    if not sample: return None
    return sample.get("ct", {}).get(ch_name)

def analyze_plant_gmo(protocol: dict, cfg: dict) -> dict:
    samples = protocol["samples"]
    pko, oko = None, None
    test_samples = []
    
    for s in samples:
        n = s["name"].lower()
        if "пко" in n: pko = s
        elif "око" in n: oko = s
        else: test_samples.append(s)
        
    if not pko or not oko:
        return {"error": "Отсутствует ПКО или ОКО. Анализ невозможен. Проверьте названия в протоколе."}

    # ПКО
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

    # ОКО
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

    # Образцы
    sample_results = []
    for s in test_samples:
        s_plant = _get_ct(s, "Hex")
        s_35s = _get_ct(s, "Rox")
        s_nos = _get_ct(s, "Fam")
        s_vpk = _get_ct(s, "Cy5")
        
        # Растение
        plant_det = "-"
        if s_plant is not None and pko_plant is not None:
            if s_plant <= pko_plant + cfg['delta_pos']: plant_det = "+"
            elif s_plant <= pko_plant + cfg['delta_low']: plant_det = "мало"
            
        # 35S / NOS
        prom_det, nos_det = "-", "-"
        if plant_det == "+":
            if s_35s and pko_35s and pko_plant and s_plant:
                if (s_35s - s_plant) <= (pko_35s - pko_plant) + cfg['delta_pos']: prom_det = "+"
            if s_nos and pko_nos and pko_plant and s_plant:
                if (s_nos - s_plant) <= (pko_nos - pko_plant) + cfg['delta_pos']: nos_det = "+"
        elif plant_det == "мало":
            if s_35s and pko_35s and s_35s <= pko_35s + cfg['delta_pos']: prom_det = "+"
            if s_nos and pko_nos and s_nos <= pko_nos + cfg['delta_pos']: nos_det = "+"

        # ВПК
        vpk_det = "+" if (s_vpk and oko_vpk and s_vpk <= oko_vpk + cfg['vpk_delta']) or (s_vpk and not oko_vpk and s_vpk <= cfg['vpk_max']) else "-"
        
        # Выводы (из листа 2)
        conclusion, action = "Результат не определён", "Проверить данные"
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
            conclusion, action = "Нет растительной ДНК", "Запрос сырья"

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
    samples = protocol["samples"]
    pko, oko, kov = None, None, None
    test_samples = []
    
    for s in samples:
        n = s["name"].lower()
        if "пко" in n: pko = s
        elif "око" in n: oko = s
        elif "ко-в" in n or "меланж" in n or "10%" in n: kov = s
        else: test_samples.append(s)

    if not pko or not oko:
        return {"error": "Отсутствует ПКО или ОКО. Проверьте названия лунок."}
        
    # Канал: FAM→Turkey, ROX→Chicken, HEX→VPK
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
        
        # Индейка (FAM)
        turkey_res = "-"
        if s_turkey is not None and s_turkey < cfg['ct_cutoff']:
            turkey_res = "В пробе содержится ДНК индейки"
        elif s_turkey is not None and s_turkey >= cfg['ct_cutoff']:
            turkey_res = "В пробе отсутствует ДНК индейки"
            
        # Курица (ROX)
        chicken_res = "Ложноотрицательный результат"
        if kov_chicken is None or kov_chicken > 29:
             chicken_res = "Ошибка калибратора (КО-В)"
        elif s_chicken is None:
             pass 
        elif s_chicken >= cfg['ct_cutoff']:
             chicken_res = "В пробе отсутствует ДНК курицы."
        elif s_chicken < kov_chicken:
             chicken_res = "В пробе содержится ДНК курицы. Количество ДНК соответствует наличию мяса курицы."
        elif s_chicken < cfg['ct_cutoff']:
             chicken_res = "В пробе содержится ДНК курицы. Количество ДНК менее ДНК 10% меланжа. Мясо курицы (продукт убоя) отсутствует."

        # ВПК
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
# ИНТЕРФЕЙС STREAMLIT (Восстановлен дизайн карточек из app1rt.py)
# ═══════════════════════════════════════════════════════════════════

st.markdown("""
<div class="header-card">
    <h1>🧬 RT-PCR Анализатор</h1>
    <p>Автоматическая интерпретация результатов из протоколов амплификатора</p>
</div>
""", unsafe_allow_html=True)

# Боковая панель: Настройки
st.sidebar.header("⚙️ Настройки порогов")
test_type_manual = st.sidebar.selectbox("Выберите тест-систему:", ["Автоопределение", "Растение/35S+FMV/NOS", "Курица/Индейка"])

st.sidebar.subheader("🌱 ГМО Скрининг")
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

# Основная область: Загрузка
uploaded = st.file_uploader("📂 Загрузите протокол амплификатора (.rtf, .txt)", type=["rtf", "txt"])

if uploaded:
    try:
        protocol = parse_rtf_protocol(uploaded.read())
    except Exception as e:
        st.error(f"Ошибка чтения файла: {e}")
        st.stop()
        
    if not protocol["samples"]:
        st.warning("В файле не найдены результаты образцов. Убедитесь, что таблица содержит колонки Fam, Hex, Rox и т.д.")
        st.stop()

    detected_name = protocol["test_name"].lower()
    if test_type_manual != "Автоопределение":
        active_test = "gmo" if "Растение" in test_type_manual else "meat"
    else:
        active_test = "gmo" if any(x in detected_name for x in ["растение", "35s", "fmv", "gmo"]) else "meat" if any(x in detected_name for x in ["курица", "индейка", "kurit"]) else "unknown"

    with st.expander("📋 Сведения о протоколе", expanded=True):
        c1, c2, c3 = st.columns(3)
        c1.write(f"**Дата:** {protocol['date']}")
        c2.write(f"**Оператор:** {protocol['operator']}")
        c3.write(f"**Тест в файле:** {protocol['test_name']}")

    if active_test == "unknown":
        st.error("Не удалось определить тип набора. Выберите его вручную в панели слева.")
        st.stop()
        
    st.divider()
    df_export = pd.DataFrame()
    
    # ==========================================
    # Отрисовка результатов ГМО
    # ==========================================
    if active_test == "gmo":
        st.markdown("### 🌱 Результаты скрининга ГМО")
        res = analyze_plant_gmo(protocol, cfg_gmo)
        if "error" in res:
            st.error(res["error"])
        else:
            # ── Контроли ──
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
                
            # ── Таблица ──
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

            # ── Детальные выводы (карточки как в app1rt.py) ──
            st.markdown("#### Подробные выводы")
            for sr in res["samples"]:
                if "ГМ растение" in sr["conclusion"]: icon, box = "🔴", "error-box"
                elif "не обнаружена (реакция валидна)" in sr["conclusion"] or "не обнаружены" in sr["conclusion"]: icon, box = "🟢", "success-box"
                elif "Нет растительной" in sr["conclusion"]: icon, box = "🟡", "warn-box"
                else: icon, box = "🔵", "info-box"
                
                st.markdown(f"""
                <div class="{box}">
                    {icon} <b>{sr['name']}</b> (лунка {sr['well']})<br>
                    {sr['conclusion']}<br>
                    <i>Действие: {sr['action']}</i>
                </div>
                """, unsafe_allow_html=True)

    # ==========================================
    # Отрисовка результатов Курица/Индейка
    # ==========================================
    elif active_test == "meat":
        st.markdown("### 🍗 Результаты определения курицы / индейки")
        res = analyze_chicken_turkey(protocol, cfg_meat)
        if "error" in res:
            st.error(res["error"])
        else:
            # ── Контроли ──
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

            # ── Таблица ──
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

            # ── Детальные выводы ──
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
                
                st.markdown(f"""
                <div class="{box}">
                    {icon} <b>{sr['name']}</b> (лунка {sr['well']})<br>
                    Индейка: {sr['turkey_result']}<br>
                    Курица: {sr['chicken_result']}
                </div>
                """, unsafe_allow_html=True)

    # Кнопка скачивания Excel (доступна для обоих наборов)
    if not df_export.empty:
        st.markdown("<br>", unsafe_allow_html=True)
        excel_data = to_excel(df_export)
        st.download_button(
            label="💾 Скачать отчет (Excel)",
            data=excel_data,
            file_name=f"Отчет_{protocol['date'].replace(' ', '_')}.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )

else:
    # Инструкция на главной странице (как в app1rt.py)
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