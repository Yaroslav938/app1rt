import streamlit as st
import re
import math
from io import BytesIO

try:
    from striprtf.striprtf import rtf_to_text
except ImportError:
    rtf_to_text = None

st.set_page_config(
    page_title="RT-PCR Анализатор",
    page_icon="🧬",
    layout="wide",
)

# ─── Custom CSS ───
st.markdown("""
<style>
    .main .block-container { max-width: 1100px; padding-top: 2rem; }
    .stTabs [data-baseweb="tab-list"] { gap: 8px; }
    .stTabs [data-baseweb="tab"] {
        padding: 10px 24px;
        font-weight: 600;
    }
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
    .info-box {
        background: #e3f2fd; border-left: 4px solid #1976d2;
        padding: 12px 16px; border-radius: 0 8px 8px 0; margin: 8px 0;
    }
    .warn-box {
        background: #fff3e0; border-left: 4px solid #f57c00;
        padding: 12px 16px; border-radius: 0 8px 8px 0; margin: 8px 0;
    }
    .error-box {
        background: #ffebee; border-left: 4px solid #d32f2f;
        padding: 12px 16px; border-radius: 0 8px 8px 0; margin: 8px 0;
    }
    .success-box {
        background: #e8f5e9; border-left: 4px solid #388e3c;
        padding: 12px 16px; border-radius: 0 8px 8px 0; margin: 8px 0;
    }
</style>
""", unsafe_allow_html=True)


# ═══════════════════════════════════════════════════════════════════
# RTF PARSER
# ═══════════════════════════════════════════════════════════════════

def parse_rtf_protocol(content_bytes: bytes) -> dict:
    """Parse an RTF protocol file from the amplifier and return structured data."""
    # Try various encodings
    text = None
    for enc in ("cp1251", "utf-8", "latin-1"):
        try:
            raw = content_bytes.decode(enc)
            break
        except UnicodeDecodeError:
            raw = None

    if raw is None:
        raise ValueError("Не удалось декодировать файл. Попробуйте другую кодировку.")

    if rtf_to_text is not None and raw.strip().startswith("{\\rtf"):
        text = rtf_to_text(raw)
    else:
        text = raw  # fallback: might be plain-text export

    result = {
        "date": "",
        "operator": "",
        "protocol_number": "",
        "result_file": "",
        "test_name": "",
        "amplification_program": "",
        "samples": [],
        "thresholds": {},
    }

    lines = text.split("\n")
    for line in lines:
        if "Дата:" in line:
            m = re.search(r'Дата:\|*\|*(.+?)(?:\||\n|$)', line)
            if m:
                result["date"] = m.group(1).strip().strip("|").strip()
        if "Оператор:" in line:
            m = re.search(r'Оператор:\|*\|*(.+?)(?:\||\n|$)', line)
            if m:
                result["operator"] = m.group(1).strip().strip("|").strip()
        if "Номер протокола:" in line:
            m = re.search(r'Номер протокола:\|*\|*(.+?)(?:\||\n|$)', line)
            if m:
                result["protocol_number"] = m.group(1).strip().strip("|").strip()
        if "Файл с результатами:" in line:
            m = re.search(r'Файл с результатами:\|*(.+?)(?:\||\n|$)', line)
            if m:
                result["result_file"] = m.group(1).strip().strip("|").strip()
        if "Тест:" in line:
            m = re.search(r'Тест:\|*\|*(.+?)(?:\||\n|$)', line)
            if m:
                result["test_name"] = m.group(1).strip().strip("|").strip()
        if "Программа амплификации:" in line:
            m = re.search(r'Программа амплификации:\|*\|*(.+?)(?:\||\n|$)', line)
            if m:
                result["amplification_program"] = m.group(1).strip().strip("|").strip()
        # Thresholds
        if "Threshold_FAM" in line:
            m = re.search(r'Threshold_FAM\s*=\s*([\d,\.]+)', line)
            if m:
                result["thresholds"]["FAM"] = float(m.group(1).replace(",", "."))
        if "Threshold_HEX" in line:
            m = re.search(r'Threshold_HEX\s*=\s*([\d,\.]+)', line)
            if m:
                result["thresholds"]["HEX"] = float(m.group(1).replace(",", "."))
        if "Threshold_ROX" in line:
            m = re.search(r'Threshold_ROX\s*=\s*([\d,\.]+)', line)
            if m:
                result["thresholds"]["ROX"] = float(m.group(1).replace(",", "."))
        if "Threshold_CY5" in line or "Threshold_Cy5" in line:
            m = re.search(r'Threshold_[Cc][Yy]5\s*=\s*([\d,\.]+)', line)
            if m:
                result["thresholds"]["Cy5"] = float(m.group(1).replace(",", "."))

    # Parse sample table
    # Find header row with channels
    channels_order = []  # e.g. ['Fam', 'Hex', 'Rox', 'Cy5', 'Cy5.5']
    sample_section = False
    for line in lines:
        if "Номер лунки" in line and "Идентификатор" in line:
            sample_section = True
            # Parse channel order from header
            parts = [p.strip() for p in line.split("|") if p.strip()]
            for p in parts:
                p_low = p.lower().replace(" ", "")
                if "fam" in p_low:
                    channels_order.append("Fam")
                elif "hex" in p_low:
                    channels_order.append("Hex")
                elif "rox" in p_low:
                    channels_order.append("Rox")
                elif "cy5.5" in p_low:
                    channels_order.append("Cy5.5")
                elif "cy5" in p_low:
                    channels_order.append("Cy5")
            continue

        if sample_section and line.strip().startswith("|") and not line.strip().startswith("|*"):
            parts = [p.strip() for p in line.split("|")]
            parts = [p for p in parts if p != ""]
            if len(parts) < 2:
                continue

            # First part is well position (A1, A2, etc.)
            well = parts[0] if re.match(r'^[A-H]\d+$', parts[0]) else ""
            # Second part is sample identifier
            sample_id_raw = parts[1] if len(parts) > 1 else ""

            # Extract sample name (remove test name in parentheses)
            sample_name = re.sub(r'\s*\(.*?\)\s*$', '', sample_id_raw).strip()

            # Remaining parts are Ct values in channel order
            ct_values = {}
            ct_parts = parts[2:] if well else parts[2:]

            # The RTF table may have empty cells between values
            # We need to match values to channels based on position
            # Re-parse the raw line to get positional data
            raw_parts = line.split("|")
            # Count cells
            # The format is: ||well|name|Fam|Hex|Rox|Cy5|Cy5.5|...
            # Find numeric values and map to channels
            numeric_vals = []
            for p in raw_parts[2:]:  # skip empty first parts
                p = p.strip()
                if re.match(r'^\d+[,\.]\d+$', p):
                    numeric_vals.append(float(p.replace(",", ".")))
                elif p == "" or p == "-":
                    numeric_vals.append(None)
                # skip non-numeric non-empty strings (like sample name)

            # Now we need a better approach: count positions after sample name
            # Let's re-parse more carefully
            ct_values = _parse_ct_from_line(line, channels_order)

            if well or sample_name:
                result["samples"].append({
                    "well": well,
                    "name": sample_name,
                    "ct": ct_values,
                })

        if sample_section and "|*" in line:
            break  # end of sample table

    return result


def _parse_ct_from_line(line: str, channels: list) -> dict:
    """Parse Ct values from a table line, mapping to channels by position."""
    ct = {}
    # Split by | but keep track of positions
    parts = line.split("|")
    
    # Find the well position and sample name first
    # Structure: ||well|name|ch1|ch2|ch3|ch4|ch5|...|
    non_empty_idx = []
    for i, p in enumerate(parts):
        if p.strip():
            non_empty_idx.append(i)

    if len(non_empty_idx) < 2:
        return ct

    # After well and name, the remaining cells map to channels
    # But we need positional mapping (empty cells = no value for that channel)
    
    # Find where the sample name ends
    # The well is always in a consistent position
    well_found = False
    name_found = False
    channel_idx = 0
    
    # Skip initial empty parts, then well, then name
    state = "pre"
    ct_cells = []
    
    for i, p in enumerate(parts):
        ps = p.strip()
        if state == "pre":
            if re.match(r'^[A-H]\d+$', ps):
                state = "well"
                continue
        elif state == "well":
            if ps and not re.match(r'^\d+[,\.]\d+$', ps):
                state = "name"
                continue
            elif re.match(r'^\d+[,\.]\d+$', ps):
                state = "ct"
                ct_cells.append(ps)
                continue
            else:
                continue
        elif state == "name":
            state = "ct"
            ct_cells.append(ps)
            continue
        elif state == "ct":
            ct_cells.append(ps)

    # Map ct_cells to channels
    for i, ch in enumerate(channels):
        if i < len(ct_cells):
            val = ct_cells[i].strip()
            if val and re.match(r'^\d+[,\.]\d+$', val):
                ct[ch] = float(val.replace(",", "."))
            else:
                ct[ch] = None
        else:
            ct[ch] = None

    return ct


# ═══════════════════════════════════════════════════════════════════
# ANALYSIS ENGINES
# ═══════════════════════════════════════════════════════════════════

def detect_test_type(protocol: dict) -> str:
    """Detect which test type based on test name in protocol."""
    name = protocol.get("test_name", "").lower()
    if "растение" in name or "35s" in name or "fmv" in name or "nos" in name:
        return "plant_gmo"
    elif "курица" in name or "индейка" in name or "kuritsa" in name:
        return "chicken_turkey"
    return "unknown"


# ─── PLANT GMO SCREENING (Растение/35S+FMV/NOS) ──────────────────

def analyze_plant_gmo(protocol: dict) -> dict:
    """
    Implements the logic from Rastenie-35S-FMV-NOS-skrining.xlt
    
    Channels mapping:
      - Yellow (Hex) → Растение (Plant DNA)
      - Orange (Rox) → П-35S+FMV (Promoter)  
      - Green (Fam)  → T-NOS (Terminator)
      - Cy5 (Red)    → ВПК (Internal amplification control)
    
    Row 9 = ПКО-С1 (Positive Control), Row 10 = ОКО (Negative Control)
    Rows 11+ = samples
    """
    samples = protocol["samples"]
    if len(samples) < 2:
        return {"error": "Недостаточно образцов. Нужен минимум ПКО и ОКО."}
    
    # Identify controls and samples
    pko = None
    oko = None
    test_samples = []
    
    for s in samples:
        name_lower = s["name"].lower()
        if "пко" in name_lower:
            pko = s
        elif "око" in name_lower:
            oko = s
        else:
            test_samples.append(s)
    
    if pko is None:
        return {"error": "Не найден ПКО (положительный контроль)."}
    if oko is None:
        return {"error": "Не найден ОКО (отрицательный контроль)."}
    
    # Extract Ct values for controls
    # Yellow/Hex → Plant, Rox/Orange → 35S+FMV, Fam/Green → NOS, Cy5/Red → VPK
    pko_plant = _get_ct(pko, "plant")    # Yellow/Hex
    pko_35s   = _get_ct(pko, "promoter") # Rox/Orange  
    pko_nos   = _get_ct(pko, "nos")      # Fam/Green
    pko_vpk   = _get_ct(pko, "vpk")      # Cy5
    
    oko_vpk   = _get_ct(oko, "vpk")      # Cy5
    
    # ─── Validate PKO ───
    # PKO is accepted if: plant Ct <= 23, promoter Ct <= 36, NOS Ct <= 36
    pko_status = "Принимается"
    pko_action = ""
    if pko_plant is None or pko_plant > 23 or \
       pko_35s is None or pko_35s > 36 or \
       pko_nos is None or pko_nos > 36:
        pko_status = "Не принимается ПКО"
        pko_action = "Повторить анализ"
    
    # ─── Validate OKO ───
    # OKO: check if any positive signals in channels + VPK
    oko_plant_ct = _get_ct(oko, "plant")
    oko_35s_ct = _get_ct(oko, "promoter")
    oko_nos_ct = _get_ct(oko, "nos")
    
    # OKO plant detection
    oko_plant_det = _detect_plant_pko(oko_plant_ct, pko_plant)
    oko_35s_det = _detect_35s_pko(oko_35s_ct, pko_35s)
    oko_nos_det = _detect_nos_pko(oko_nos_ct, pko_nos)
    
    oko_has_signal = oko_plant_det in ("+", "мало") or oko_35s_det == "+" or oko_nos_det == "+"
    oko_vpk_ok = oko_vpk is not None and oko_vpk < (pko_vpk + 5 if pko_vpk else 36)
    
    oko_status = "Принимается"
    oko_action = ""
    if oko_has_signal and oko_vpk_ok:
        # Signal number = 9+10=19 → Контаминация
        oko_status = "Контаминация"
        oko_action = "Провести деконтаминацию"
    elif oko_has_signal and not oko_vpk_ok:
        oko_status = "Контаминация"
        oko_action = "Провести деконтаминацию"
    elif not oko_has_signal and not oko_vpk_ok:
        oko_status = "Не принимается ВПК"
        oko_action = "Повторить анализ"
    # else: accepted
    
    # ─── Analyze samples ───
    sample_results = []
    for s in test_samples:
        s_plant = _get_ct(s, "plant")
        s_35s   = _get_ct(s, "promoter")
        s_nos   = _get_ct(s, "nos")
        s_vpk   = _get_ct(s, "vpk")
        
        # Plant detection (Yellow/Hex)
        # PKO row 9: +3.4 threshold, +10.5 for "мало"
        plant_det = _detect_plant(s_plant, pko_plant)
        
        # 35S+FMV detection (Orange/Rox) - depends on plant detection
        promoter_det = _detect_promoter_sample(plant_det, s_35s, s_plant, pko_35s, pko_plant)
        
        # NOS detection (Green/Fam) - depends on plant detection 
        nos_det = _detect_nos_sample(plant_det, s_nos, s_plant, pko_nos, pko_plant)
        
        # VPK detection (Cy5/Red)
        vpk_det = _detect_vpk(s_vpk, oko_vpk)
        
        # Final conclusion
        conclusion, action = _plant_conclusion(plant_det, promoter_det, nos_det, vpk_det)
        
        sample_results.append({
            "name": s["name"],
            "well": s["well"],
            "plant_ct": s_plant,
            "plant_det": plant_det,
            "promoter_ct": s_35s,
            "promoter_det": promoter_det,
            "nos_ct": s_nos,
            "nos_det": nos_det,
            "vpk_ct": s_vpk,
            "vpk_det": vpk_det,
            "conclusion": conclusion,
            "action": action,
        })
    
    return {
        "type": "plant_gmo",
        "pko": {
            "name": pko["name"],
            "well": pko["well"],
            "plant_ct": pko_plant,
            "promoter_ct": pko_35s,
            "nos_ct": pko_nos,
            "vpk_ct": pko_vpk,
            "status": pko_status,
            "action": pko_action,
        },
        "oko": {
            "name": oko["name"],
            "well": oko["well"],
            "vpk_ct": oko_vpk,
            "status": oko_status,
            "action": oko_action,
        },
        "samples": sample_results,
    }


def _get_ct(sample: dict, target: str) -> float | None:
    """Get Ct value for a target from a sample's channels."""
    ct = sample.get("ct", {})
    if target == "plant":    # Yellow / Hex
        return ct.get("Hex") or ct.get("Yellow")
    elif target == "promoter":  # Orange / Rox
        return ct.get("Rox") or ct.get("Orange")
    elif target == "nos":       # Green / Fam
        return ct.get("Fam") or ct.get("Green")
    elif target == "vpk":       # Red / Cy5
        return ct.get("Cy5") or ct.get("Red")
    return None


def _detect_plant_pko(ct, pko_ct):
    """Plant detection for PKO row."""
    if ct is None:
        return "-"
    if pko_ct is None:
        return "+"  if ct > 0 else "-"
    if ct > pko_ct + 3.4:
        if ct > pko_ct + 10.5:
            return "-"
        return "мало"
    return "+"

def _detect_35s_pko(ct, pko_ct):
    """35S detection for PKO/OKO."""
    if ct is None:
        return "-"
    if pko_ct is None:
        return "+"
    return "+" if ct <= pko_ct + 3.4 else "-"

def _detect_nos_pko(ct, pko_ct):
    """NOS detection for PKO/OKO."""
    if ct is None:
        return "-"
    if pko_ct is None:
        return "+"
    return "+" if ct <= pko_ct + 3.4 else "-"


def _detect_plant(ct, pko_ct):
    """Plant detection for samples: compare to PKO plant Ct."""
    if ct is None:
        return "-"
    if pko_ct is None:
        return "-"
    if ct > pko_ct + 3.4:
        if ct > pko_ct + 10.5:
            return "-"
        return "мало"
    return "+"


def _detect_promoter_sample(plant_det, s_35s, s_plant, pko_35s, pko_plant):
    """35S+FMV detection for samples, depends on plant detection result."""
    if plant_det == "+":
        if s_35s is None:
            return "-"
        if pko_35s is not None and pko_plant is not None and s_plant is not None:
            delta_threshold = (pko_35s - pko_plant) + 3.4
            if (s_35s - s_plant) > delta_threshold:
                return "-"
            return "+"
        return "-"
    elif plant_det == "мало":
        if s_35s is None:
            return "-"
        if pko_35s is not None:
            return "+" if s_35s <= pko_35s + 3.4 else "-"
        return "-"
    return "-"


def _detect_nos_sample(plant_det, s_nos, s_plant, pko_nos, pko_plant):
    """NOS detection for samples, depends on plant detection result."""
    if plant_det == "+":
        if s_nos is None:
            return "-"
        if pko_nos is not None and pko_plant is not None and s_plant is not None:
            delta_threshold = (pko_nos - pko_plant) + 3.4
            if (s_nos - s_plant) > delta_threshold:
                return "-"
            return "+"
        return "-"
    elif plant_det == "мало":
        if s_nos is None:
            return "-"
        if pko_nos is not None:
            return "+" if s_nos <= pko_nos + 3.4 else "-"
        return "-"
    return "-"


def _detect_vpk(ct, oko_vpk):
    """VPK (internal control) detection."""
    if ct is None:
        return "-"
    threshold = (oko_vpk + 5) if oko_vpk is not None else 36
    return "+" if ct <= threshold else "-"


def _plant_conclusion(plant_det, promoter_det, nos_det, vpk_det):
    """
    Final conclusion logic from Лист2:
    Code 1: plant+ AND (promoter+ OR nos+) → GM plant, high DNA
    Code 2: plant+ AND promoter- AND nos- → No promoters/terminators found
    Code 3: plant="мало" AND promoter- AND nos- → No promoters/terminators found
    Code 4: plant- AND vpk+ → False negative
    Code 5: plant- AND all negative AND vpk- → No plant DNA
    Code 6: plant="мало" AND (promoter+ OR nos+) → GM plant, low DNA
    """
    if plant_det == "+" and (promoter_det == "+" or nos_det == "+"):
        return ("ГМ растение; количество ДНК выше предела количественного определения",
                "Идентификация и определение количества %")
    
    if plant_det == "+" and promoter_det == "-" and nos_det == "-":
        return ("Промоторы 35S, FMV и терминатор NOS не обнаружены",
                "Обнаружение сои линий BPS-CV-127-9 и MON87701")
    
    if plant_det == "мало" and promoter_det == "-" and nos_det == "-":
        return ("Промоторы 35S, FMV и терминатор NOS не обнаружены",
                "Концентрирование ДНК; запрос сырья")
    
    if plant_det == "-" and vpk_det == "+":
        return ("ДНК растения не обнаружена (реакция валидна)",
                "")
    
    if plant_det == "мало" and (promoter_det == "+" or nos_det == "+"):
        return ("ГМ растение; для точного количественного анализа требуется концентрирование ДНК",
                "Концентрирование и идентификация")
    
    if plant_det == "-" and promoter_det == "-" and nos_det == "-" and vpk_det == "-":
        return ("Нет растительной ДНК",
                "Запрос сырья")
    
    return ("Результат не определён", "Проверить данные")


# ─── CHICKEN / TURKEY DETECTION (Курица-Индейка) ─────────────────

def analyze_chicken_turkey(protocol: dict) -> dict:
    """
    Implements the logic from Kuritsa-Indeika.xlt
    
    Channels:
      - FAM/Green → Turkey DNA
      - ROX/Orange → Chicken DNA
      - HEX/Yellow → VPK (Internal control)
    
    Row 9 = ПКО (Positive Control)
    Row 10 = КО-В меланж 10% (Calibration sample)
    Row 11 = ОКО (Negative Control)
    Rows 12+ = samples
    """
    samples = protocol["samples"]
    
    pko = None
    kov = None  # КО-В меланж 10%
    oko = None
    test_samples = []
    
    for s in samples:
        name_lower = s["name"].lower()
        if "пко" in name_lower:
            pko = s
        elif "ко-в" in name_lower or "меланж" in name_lower:
            kov = s
        elif "око" in name_lower:
            oko = s
        else:
            test_samples.append(s)
    
    if pko is None:
        return {"error": "Не найден ПКО (положительный контроль)."}
    if oko is None:
        return {"error": "Не найден ОКО (отрицательный контроль)."}
    
    # Channel mapping for this test: FAM→Turkey, ROX→Chicken, HEX→VPK
    pko_turkey  = _get_ct_chicken(pko, "turkey")
    pko_chicken = _get_ct_chicken(pko, "chicken")
    pko_vpk     = _get_ct_chicken(pko, "vpk")
    
    # PKO validation
    pko_turkey_det = _det_35(pko_turkey)
    pko_chicken_det = _det_35(pko_chicken)
    pko_vpk_det = _det_35(pko_vpk)
    
    # KO-V (calibration)
    kov_results = None
    if kov:
        kov_chicken = _get_ct_chicken(kov, "chicken")
        kov_vpk = _get_ct_chicken(kov, "vpk")
        # KO-V chicken: check if > 29 → "-", else "+"
        kov_chicken_det = "-"
        if kov_chicken is not None:
            kov_chicken_det = "-" if kov_chicken > 29 else "+"
        kov_vpk_det = _det_35(kov_vpk)
        kov_results = {
            "name": kov["name"],
            "well": kov["well"],
            "chicken_ct": kov_chicken,
            "chicken_det": kov_chicken_det,
            "vpk_ct": kov_vpk,
            "vpk_det": kov_vpk_det,
        }
    
    # OKO
    oko_turkey  = _get_ct_chicken(oko, "turkey")
    oko_chicken = _get_ct_chicken(oko, "chicken")
    oko_vpk     = _get_ct_chicken(oko, "vpk")
    oko_turkey_det = _det_35(oko_turkey)
    oko_chicken_det = _det_35(oko_chicken)
    oko_vpk_det = _det_35(oko_vpk)
    
    # Analyze samples
    # The KO-V calibration chicken Ct is reference for "количество соответствует наличию мяса"
    kov_chicken_ct = _get_ct_chicken(kov, "chicken") if kov else None
    kov_chicken_det_val = kov_results["chicken_det"] if kov_results else "-"
    
    sample_results = []
    for s in test_samples:
        s_turkey  = _get_ct_chicken(s, "turkey")
        s_chicken = _get_ct_chicken(s, "chicken")
        s_vpk     = _get_ct_chicken(s, "vpk")
        
        # Turkey detection: Ct < 35 → detected
        turkey_result = "-"
        if s_turkey is not None and s_turkey < 35:
            turkey_result = "В пробе содержится ДНК индейки"
        elif s_turkey is not None and s_turkey >= 35:
            turkey_result = "В пробе отсутствует ДНК индейки"
        else:
            turkey_result = "-"
        
        # Chicken detection: depends on KO-V calibration
        chicken_result = _chicken_conclusion(s_chicken, kov_chicken_ct, kov_chicken_det_val)
        
        # VPK
        vpk_det = _det_35(s_vpk)
        
        sample_results.append({
            "name": s["name"],
            "well": s["well"],
            "turkey_ct": s_turkey,
            "turkey_result": turkey_result,
            "chicken_ct": s_chicken,
            "chicken_result": chicken_result,
            "vpk_ct": s_vpk,
            "vpk_det": vpk_det,
        })
    
    return {
        "type": "chicken_turkey",
        "pko": {
            "name": pko["name"],
            "well": pko["well"],
            "turkey_ct": pko_turkey,
            "turkey_det": pko_turkey_det,
            "chicken_ct": pko_chicken,
            "chicken_det": pko_chicken_det,
            "vpk_ct": pko_vpk,
            "vpk_det": pko_vpk_det,
        },
        "kov": kov_results,
        "oko": {
            "name": oko["name"],
            "well": oko["well"],
            "turkey_ct": oko_turkey,
            "turkey_det": oko_turkey_det,
            "chicken_ct": oko_chicken,
            "chicken_det": oko_chicken_det,
            "vpk_ct": oko_vpk,
            "vpk_det": oko_vpk_det,
        },
        "samples": sample_results,
    }


def _get_ct_chicken(sample: dict, target: str) -> float | None:
    """Channel mapping for Chicken-Turkey test."""
    ct = sample.get("ct", {})
    if target == "turkey":   # FAM / Green
        return ct.get("Fam") or ct.get("Green")
    elif target == "chicken": # ROX / Orange
        return ct.get("Rox") or ct.get("Orange")
    elif target == "vpk":     # HEX / Yellow
        return ct.get("Hex") or ct.get("Yellow")
    return None


def _det_35(ct):
    """Simple detection: Ct < 35 → '+', else '-'."""
    if ct is not None and ct < 35:
        return "+"
    return "-"


def _chicken_conclusion(s_chicken, kov_ct, kov_det):
    """Chicken detection conclusion per Excel logic."""
    if kov_det != "+":
        return "Ложноотрицательный результат"
    
    if s_chicken is None:
        return "-"
    
    if s_chicken >= 35:
        return "В пробе отсутствует ДНК курицы."
    
    if kov_ct is not None and s_chicken < kov_ct:
        return "В пробе содержится ДНК курицы. Количество ДНК соответствует наличию мяса курицы."
    
    if s_chicken < 35:
        return "В пробе содержится ДНК курицы. Количество ДНК менее ДНК 10% меланжа. Мясо курицы (продукт убоя) отсутствует."
    
    return "В пробе отсутствует ДНК курицы."


# ═══════════════════════════════════════════════════════════════════
# STREAMLIT UI
# ═══════════════════════════════════════════════════════════════════

st.markdown("""
<div class="header-card">
    <h1>🧬 RT-PCR Анализатор</h1>
    <p>Автоматическая интерпретация результатов real-time ПЦР из протоколов амплификатора</p>
</div>
""", unsafe_allow_html=True)

uploaded = st.file_uploader(
    "Загрузите файл протокола (.rtf или .txt)",
    type=["rtf", "txt", "doc"],
    help="Файл протокола из амплификатора (формат RTF или текстовый)"
)

if uploaded is not None:
    content = uploaded.read()
    
    try:
        protocol = parse_rtf_protocol(content)
    except Exception as e:
        st.error(f"Ошибка при чтении файла: {e}")
        st.stop()
    
    if not protocol["samples"]:
        st.error("Не удалось извлечь данные из файла. Проверьте формат протокола.")
        st.stop()
    
    # ─── Protocol metadata ───
    with st.expander("📋 Информация о протоколе", expanded=True):
        col1, col2 = st.columns(2)
        with col1:
            st.markdown(f"**Дата:** {protocol['date']}")
            st.markdown(f"**Оператор:** {protocol['operator']}")
            st.markdown(f"**Тест:** {protocol['test_name']}")
        with col2:
            st.markdown(f"**Протокол №:** {protocol['protocol_number']}")
            st.markdown(f"**Файл:** {protocol['result_file']}")
            st.markdown(f"**Программа:** {protocol['amplification_program']}")
    
    # ─── Raw data preview ───
    with st.expander("📊 Исходные данные Ct", expanded=False):
        import pandas as pd
        raw_rows = []
        for s in protocol["samples"]:
            row = {"Лунка": s["well"], "Образец": s["name"]}
            for ch, val in s["ct"].items():
                row[f"Ct ({ch})"] = val if val is not None else "—"
            raw_rows.append(row)
        if raw_rows:
            st.dataframe(pd.DataFrame(raw_rows), use_container_width=True, hide_index=True)
    
    # ─── Detect test type and analyze ───
    test_type = detect_test_type(protocol)
    
    # Allow manual override
    type_options = {
        "plant_gmo": "Растение/35S+FMV/NOS скрининг",
        "chicken_turkey": "Курица/Индейка",
    }
    
    if test_type == "unknown":
        st.warning("Тип теста не определён автоматически. Выберите вручную:")
        test_type = st.selectbox(
            "Тип анализа",
            options=list(type_options.keys()),
            format_func=lambda x: type_options[x],
        )
    else:
        st.markdown(f'<div class="info-box">Определён тип анализа: <b>{type_options.get(test_type, test_type)}</b></div>', unsafe_allow_html=True)
    
    st.divider()
    
    # ═══ PLANT GMO RESULTS ═══
    if test_type == "plant_gmo":
        results = analyze_plant_gmo(protocol)
        
        if "error" in results:
            st.error(results["error"])
            st.stop()
        
        st.markdown("### 🌱 Результаты скрининга ГМО")
        
        # Controls status
        st.markdown("#### Контроли")
        c1, c2 = st.columns(2)
        
        with c1:
            pko = results["pko"]
            box_class = "success-box" if pko["status"] == "Принимается" else "error-box"
            st.markdown(f"""
            <div class="{box_class}">
                <b>ПКО ({pko['name']})</b> — лунка {pko['well']}<br>
                Растение: <code>{pko['plant_ct'] or '—'}</code> | 
                35S+FMV: <code>{pko['promoter_ct'] or '—'}</code> | 
                NOS: <code>{pko['nos_ct'] or '—'}</code> | 
                ВПК: <code>{pko['vpk_ct'] or '—'}</code><br>
                <b>Статус: {pko['status']}</b> {(' — ' + pko['action']) if pko['action'] else ''}
            </div>
            """, unsafe_allow_html=True)
        
        with c2:
            oko = results["oko"]
            box_class = "success-box" if oko["status"] == "Принимается" else "error-box"
            st.markdown(f"""
            <div class="{box_class}">
                <b>ОКО ({oko['name']})</b> — лунка {oko['well']}<br>
                ВПК: <code>{oko['vpk_ct'] or '—'}</code><br>
                <b>Статус: {oko['status']}</b> {(' — ' + oko['action']) if oko['action'] else ''}
            </div>
            """, unsafe_allow_html=True)
        
        # Sample results
        st.markdown("#### Образцы")
        
        import pandas as pd
        rows = []
        for sr in results["samples"]:
            rows.append({
                "Лунка": sr["well"],
                "Образец": sr["name"],
                "Растение (Ct)": f"{sr['plant_ct']:.1f}" if sr['plant_ct'] else "—",
                "Растение": sr["plant_det"],
                "35S+FMV (Ct)": f"{sr['promoter_ct']:.1f}" if sr['promoter_ct'] else "—",
                "35S+FMV": sr["promoter_det"],
                "NOS (Ct)": f"{sr['nos_ct']:.1f}" if sr['nos_ct'] else "—",
                "NOS": sr["nos_det"],
                "ВПК (Ct)": f"{sr['vpk_ct']:.1f}" if sr['vpk_ct'] else "—",
                "ВПК": sr["vpk_det"],
                "Вывод": sr["conclusion"],
                "Действие": sr["action"],
            })
        
        if rows:
            df = pd.DataFrame(rows)
            
            def highlight_detection(val):
                if val == "+":
                    return "background-color: #ffcdd2; font-weight: bold"
                elif val == "мало":
                    return "background-color: #fff9c4; font-weight: bold"
                elif val == "-":
                    return "background-color: #c8e6c9"
                return ""
            
            det_cols = ["Растение", "35S+FMV", "NOS", "ВПК"]
            styled = df.style.applymap(highlight_detection, subset=det_cols)
            st.dataframe(styled, use_container_width=True, hide_index=True, height=min(400, 40 + 35 * len(rows)))
        
        # Detailed conclusions
        st.markdown("#### Подробные выводы")
        for sr in results["samples"]:
            if "ГМ растение" in sr["conclusion"]:
                icon = "🔴"
                box = "error-box"
            elif "не обнаружена (реакция валидна)" in sr["conclusion"]:
                icon = "🟢"
                box = "success-box"
            elif "Промоторы" in sr["conclusion"] and "не обнаружены" in sr["conclusion"]:
                icon = "🟢"
                box = "success-box"
            elif "Нет растительной" in sr["conclusion"]:
                icon = "🟡"
                box = "warn-box"
            else:
                icon = "🔵"
                box = "info-box"
            
            st.markdown(f"""
            <div class="{box}">
                {icon} <b>{sr['name']}</b> (лунка {sr['well']})<br>
                {sr['conclusion']}<br>
                <i>Действие: {sr['action']}</i>
            </div>
            """, unsafe_allow_html=True)
    
    # ═══ CHICKEN/TURKEY RESULTS ═══
    elif test_type == "chicken_turkey":
        results = analyze_chicken_turkey(protocol)
        
        if "error" in results:
            st.error(results["error"])
            st.stop()
        
        st.markdown("### 🍗 Результаты определения курицы / индейки")
        
        # Controls
        st.markdown("#### Контроли")
        cols = st.columns(3 if results.get("kov") else 2)
        
        with cols[0]:
            pko = results["pko"]
            st.markdown(f"""
            <div class="info-box">
                <b>ПКО ({pko['name']})</b> — {pko['well']}<br>
                Индейка (FAM): <code>{pko['turkey_ct'] or '—'}</code> [{pko['turkey_det']}]<br>
                Курица (ROX): <code>{pko['chicken_ct'] or '—'}</code> [{pko['chicken_det']}]<br>
                ВПК (HEX): <code>{pko['vpk_ct'] or '—'}</code> [{pko['vpk_det']}]
            </div>
            """, unsafe_allow_html=True)
        
        if results.get("kov"):
            with cols[1]:
                kov = results["kov"]
                st.markdown(f"""
                <div class="info-box">
                    <b>КО-В ({kov['name']})</b> — {kov['well']}<br>
                    Курица (ROX): <code>{kov['chicken_ct'] or '—'}</code> [{kov['chicken_det']}]<br>
                    ВПК (HEX): <code>{kov['vpk_ct'] or '—'}</code> [{kov['vpk_det']}]
                </div>
                """, unsafe_allow_html=True)
        
        with cols[-1]:
            oko = results["oko"]
            st.markdown(f"""
            <div class="info-box">
                <b>ОКО ({oko['name']})</b> — {oko['well']}<br>
                Индейка: <code>{oko['turkey_ct'] or '—'}</code> [{oko['turkey_det']}]<br>
                Курица: <code>{oko['chicken_ct'] or '—'}</code> [{oko['chicken_det']}]<br>
                ВПК: <code>{oko['vpk_ct'] or '—'}</code> [{oko['vpk_det']}]
            </div>
            """, unsafe_allow_html=True)
        
        # Samples
        st.markdown("#### Образцы")
        import pandas as pd
        rows = []
        for sr in results["samples"]:
            rows.append({
                "Лунка": sr["well"],
                "Образец": sr["name"],
                "Индейка Ct (FAM)": f"{sr['turkey_ct']:.1f}" if sr['turkey_ct'] else "—",
                "Индейка": sr["turkey_result"],
                "Курица Ct (ROX)": f"{sr['chicken_ct']:.1f}" if sr['chicken_ct'] else "—",
                "Курица": sr["chicken_result"],
                "ВПК Ct (HEX)": f"{sr['vpk_ct']:.1f}" if sr['vpk_ct'] else "—",
                "ВПК": sr["vpk_det"],
            })
        
        if rows:
            df = pd.DataFrame(rows)
            st.dataframe(df, use_container_width=True, hide_index=True, height=min(500, 40 + 35 * len(rows)))
        
        # Detailed conclusions
        st.markdown("#### Подробные выводы")
        for sr in results["samples"]:
            turkey_found = "содержится ДНК индейки" in sr["turkey_result"]
            chicken_found = "содержится ДНК курицы" in sr["chicken_result"]
            
            if turkey_found or chicken_found:
                box = "warn-box"
                icon = "🔴" if "мяса курицы" in sr["chicken_result"] else "🟡"
            elif "Ложноотрицательный" in sr["chicken_result"]:
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

else:
    # Landing page
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
        3. Приложение автоматически определит тип теста
        4. Результаты будут рассчитаны по алгоритмам тест-систем
        """)
