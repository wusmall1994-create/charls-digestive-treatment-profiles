from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
from scipy.stats import chi2, norm


REPO = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("CHARLS_DATA_DIR", REPO / "data"))
OUT = REPO / "source_data"
OUT.mkdir(parents=True, exist_ok=True)

WAVES = {1: 2011, 2: 2013, 3: 2015, 4: 2018}
INTERVALS = {(1, 2): "2011-2013", (2, 3): "2013-2015", (3, 4): "2015-2018"}
DISEASES = [
    "hibpe", "diabe", "cancre", "lunge", "hearte", "stroke", "psyche",
    "arthre", "dyslipe", "livere", "kidneye", "asthmae", "memrye",
]


def numeric(x: pd.Series) -> pd.Series:
    return pd.to_numeric(x, errors="coerce")


def normalize_id(x: pd.Series) -> pd.Series:
    s = x.astype(str).str.strip()
    # Wave 1 public files use an 11-character ID; Harmonized CHARLS inserts a
    # zero before the final two-character person identifier.
    return s.where(s.str.len().ne(11), s.str.slice(0, 9) + "0" + s.str.slice(9))


def wmean(x: pd.Series, w: pd.Series) -> float:
    ok = x.notna() & w.notna() & (w > 0)
    return float(np.average(x.loc[ok], weights=w.loc[ok])) if ok.any() else np.nan


def weighted_table(df: pd.DataFrame, row: str, col: str, weight: str) -> pd.DataFrame:
    z = df.dropna(subset=[row, col, weight]).copy()
    tab = z.pivot_table(index=row, columns=col, values=weight, aggfunc="sum", fill_value=0)
    return tab.div(tab.sum(axis=1), axis=0) * 100


def build_profiles() -> dict[int, pd.DataFrame]:
    specs = {
        1: ("w1_health.dta", ["da010_10_s1", "da010_10_s2", "da010_10_s3", "da010_10_s4"]),
        2: ("w2_health.dta", ["da010_10_s1", "da010_10_s2", "da010_10_s3", "da010_10_s4"]),
        3: ("w3_health.dta", ["da010_10_s1", "da010_10_s2", "da010_10_s3", "da010_10_s4"]),
        4: ("w4_health.dta", ["da010_w4_10__s1", "da010_w4_10__s2", "da010_w4_10__s3", "da010_w4_10__s4"]),
    }
    result: dict[int, pd.DataFrame] = {}
    for wave, (fname, cols) in specs.items():
        d = pd.read_stata(DATA / fname, columns=["ID", *cols], convert_categoricals=False)
        d["ID"] = normalize_id(d["ID"])
        selected = pd.DataFrame(index=d.index)
        for i, c in enumerate(cols, start=1):
            selected[i] = numeric(d[c]).eq(i)
        observed = d[cols].notna().any(axis=1)
        contradiction = selected[4] & selected[[1, 2, 3]].any(axis=1)
        profile = pd.Series(pd.NA, index=d.index, dtype="string")
        profile.loc[observed & ~contradiction & selected[4]] = "None listed"
        profile.loc[observed & ~contradiction & selected[1] & ~selected[2]] = "TCM without Western"
        profile.loc[observed & ~contradiction & selected[2] & ~selected[1]] = "Western without TCM"
        profile.loc[observed & ~contradiction & selected[1] & selected[2]] = "Combined TCM and Western"
        profile.loc[observed & ~contradiction & selected[3] & ~selected[1] & ~selected[2]] = "Other only"
        result[wave] = pd.DataFrame(
            {
                "ID": d["ID"],
                f"profile{wave}": profile,
                f"profile_observed{wave}": observed,
                f"profile_contradiction{wave}": contradiction,
            }
        ).drop_duplicates("ID")
    return result


def load_harmonized() -> pd.DataFrame:
    cols = ["ID", "ragender", "raeduc_c"]
    for w in WAVES:
        cols.extend(
            [
                f"r{w}iwstat", f"r{w}agey", f"r{w}digeste", f"r{w}digestf",
                f"r{w}rural2", f"h{w}rural", f"r{w}mstat", f"r{w}higov",
                f"r{w}adla_c", f"r{w}wtrespb", f"r{w}rxdigest",
                f"r{w}rxdigest_c",
            ]
        )
        cols.extend(f"r{w}{d}" for d in DISEASES)
    h = pd.read_stata(DATA / "harmonized.dta", columns=cols, convert_categoricals=False)
    h["ID"] = normalize_id(h["ID"])
    return h


def build_person_wave(h: pd.DataFrame, profiles: dict[int, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for w, year in WAVES.items():
        z = h[["ID", "ragender", "raeduc_c"] + [c for c in h if c.startswith(f"r{w}") or c == f"h{w}rural"]].copy()
        rename = {
            f"r{w}iwstat": "iwstat", f"r{w}agey": "age", f"r{w}digeste": "digest",
            f"r{w}digestf": "dispute", f"r{w}rural2": "rural_hukou",
            f"h{w}rural": "rural_residence", f"r{w}mstat": "mstat",
            f"r{w}higov": "public_insurance", f"r{w}adla_c": "adl_count",
            f"r{w}wtrespb": "weight", f"r{w}rxdigest": "western_med",
            f"r{w}rxdigest_c": "any_med",
        }
        rename.update({f"r{w}{d}": d for d in DISEASES})
        z = z.rename(columns=rename).merge(profiles[w], on="ID", how="left")
        z = z.rename(columns={f"profile{w}": "profile", f"profile_observed{w}": "profile_observed", f"profile_contradiction{w}": "profile_contradiction"})
        z["wave"] = w
        z["year"] = year
        z["female"] = numeric(z["ragender"]).eq(2).astype(float)
        edu = numeric(z["raeduc_c"])
        z["education"] = pd.Categorical(
            np.select([edu.between(1, 5), edu.between(6, 7), edu.between(8, 10)], ["Lower", "Upper secondary/vocational", "Tertiary"], default=None),
            categories=["Lower", "Upper secondary/vocational", "Tertiary"],
        )
        z["partnered"] = numeric(z["mstat"]).isin([1, 3]).astype(float)
        z.loc[numeric(z["mstat"]).isna(), "partnered"] = np.nan
        z["public_insurance"] = numeric(z["public_insurance"]).eq(1).astype(float)
        z.loc[numeric(z["public_insurance"]).isna(), "public_insurance"] = np.nan
        z["adl_difficulty"] = numeric(z["adl_count"]).gt(0).astype(float)
        z.loc[numeric(z["adl_count"]).isna(), "adl_difficulty"] = np.nan
        vals = pd.concat([numeric(z[d]).rename(d) for d in DISEASES], axis=1)
        observed_n = vals.notna().sum(axis=1)
        z["comorbidity_count"] = vals.eq(1).sum(axis=1).where(observed_n >= 10)
        z["digest"] = numeric(z["digest"])
        z["dispute"] = numeric(z["dispute"])
        z["rural_hukou"] = numeric(z["rural_hukou"])
        z["rural_residence"] = numeric(z["rural_residence"])
        z["weight"] = numeric(z["weight"])
        z["age"] = numeric(z["age"])
        z["iwstat"] = numeric(z["iwstat"])
        z["none_listed"] = z["profile"].eq("None listed")
        z["any_listed"] = z["profile"].notna() & ~z["none_listed"]
        z["joint_group"] = pd.Categorical(
            np.select(
                [
                    z["rural_hukou"].eq(0) & z["rural_residence"].eq(0),
                    z["rural_hukou"].eq(1) & z["rural_residence"].eq(0),
                    z["rural_hukou"].eq(0) & z["rural_residence"].eq(1),
                    z["rural_hukou"].eq(1) & z["rural_residence"].eq(1),
                ],
                ["Urban hukou/urban residence", "Rural hukou/urban residence", "Urban hukou/rural residence", "Rural hukou/rural residence"],
                default=None,
            ),
            categories=["Urban hukou/urban residence", "Rural hukou/urban residence", "Urban hukou/rural residence", "Rural hukou/rural residence"],
        )
        rows.append(z)
    return pd.concat(rows, ignore_index=True)


def build_intervals(pw: pd.DataFrame) -> pd.DataFrame:
    pieces = []
    keep_base = [
        "ID", "year", "age", "female", "education", "partnered", "public_insurance",
        "adl_difficulty", "comorbidity_count", "rural_hukou", "rural_residence",
        "joint_group", "weight", "digest", "dispute", "profile", "none_listed",
        "any_listed", "iwstat",
    ]
    for (w0, w1), label in INTERVALS.items():
        a = pw.loc[pw.wave.eq(w0), keep_base].copy().add_suffix("_base")
        a = a.rename(columns={"ID_base": "ID"})
        b = pw.loc[pw.wave.eq(w1), ["ID", "year", "digest", "dispute", "profile", "none_listed", "any_listed", "iwstat"]].copy().add_suffix("_follow")
        b = b.rename(columns={"ID_follow": "ID"})
        z = a.merge(b, on="ID", how="left")
        z["interval"] = label
        z["baseline_wave"] = w0
        z["follow_wave"] = w1
        z["eligible_base"] = z["digest_base"].eq(1) & z["profile_base"].notna() & z["age_base"].ge(45) & z["weight_base"].gt(0)
        z["observed_follow_profile"] = z["digest_follow"].eq(1) & z["profile_follow"].notna()
        z["known_death"] = z["iwstat_follow"].eq(5)
        z["alive_nonresponse"] = z["iwstat_follow"].eq(4)
        z["unknown_vital"] = z["iwstat_follow"].eq(9)
        z["responded_no_eligible_profile"] = z["iwstat_follow"].eq(1) & ~z["observed_follow_profile"]
        z["transition4"] = pd.Series(pd.NA, index=z.index, dtype="string")
        obs = z["observed_follow_profile"] & z["profile_base"].notna()
        z.loc[obs & z["none_listed_base"] & z["none_listed_follow"], "transition4"] = "Stable none-listed reporting"
        z.loc[obs & z["none_listed_base"] & z["any_listed_follow"], "transition4"] = "Initiation of listed-treatment reporting"
        z.loc[obs & z["any_listed_base"] & z["none_listed_follow"], "transition4"] = "Cessation of listed-treatment reporting"
        z.loc[obs & z["any_listed_base"] & z["any_listed_follow"], "transition4"] = "Continued listed-treatment reporting"
        pieces.append(z)
    return pd.concat(pieces, ignore_index=True)


def fit_gee(df: pd.DataFrame, outcome: str, weight_col: str, label: str, interaction: bool = False):
    exposure = "C(joint_group_base, Treatment(reference='Urban hukou/urban residence'))"
    interval_term = "C(interval, Treatment(reference='2011-2013'))"
    formula = (
        f"{outcome} ~ {exposure}"
        f" + {interval_term}" + (f" + {exposure}:{interval_term}" if interaction else "") + " + age_base + female_base"
        " + C(education_base, Treatment(reference='Lower')) + partnered_base"
        " + public_insurance_base + adl_difficulty_base + comorbidity_count_base"
    )
    cols = [
        outcome, "ID", "joint_group_base", "interval", "age_base", "female_base",
        "education_base", "partnered_base", "public_insurance_base", "adl_difficulty_base",
        "comorbidity_count_base", weight_col,
    ]
    z = df[cols].dropna().copy()
    z[outcome] = z[outcome].astype(int)
    z[weight_col] = z.groupby("interval", observed=True)[weight_col].transform(lambda x: x / x.mean())
    y, X = patsy.dmatrices(formula, z, return_type="dataframe")
    model = sm.GEE(
        y.iloc[:, 0], X, groups=z["ID"].to_numpy(), family=sm.families.Poisson(),
        cov_struct=sm.cov_struct.Independence(), weights=z[weight_col].to_numpy(),
    )
    result = model.fit(maxiter=200)
    params = result.params.to_numpy()
    cov = result.cov_params().to_numpy()
    design_info = X.design_info
    groups = list(z["joint_group_base"].cat.categories)

    def standardized(group: str):
        nd = z.copy()
        nd["joint_group_base"] = pd.Categorical([group] * len(nd), categories=groups)
        Xm = np.asarray(patsy.build_design_matrices([design_info], nd)[0])
        mu = np.exp(Xm @ params)
        p = float(mu.mean())
        grad = (mu[:, None] * Xm).mean(axis=0)
        se = float(math.sqrt(max(0.0, grad @ cov @ grad)))
        return p, grad, se

    def contrast_rows(eval_data: pd.DataFrame, interval_label: str):
        local = []
        def std_local(group: str):
            nd = eval_data.copy()
            nd["joint_group_base"] = pd.Categorical([group] * len(nd), categories=groups)
            Xm = np.asarray(patsy.build_design_matrices([design_info], nd)[0])
            mu = np.exp(Xm @ params)
            p = float(mu.mean())
            grad = (mu[:, None] * Xm).mean(axis=0)
            se = float(math.sqrt(max(0.0, grad @ cov @ grad)))
            return p, grad, se
        base_p, base_g, _ = std_local(groups[0])
        for group in groups:
            p, g, se = std_local(group)
            rd = p - base_p
            grd = g - base_g
            rd_se = float(math.sqrt(max(0.0, grd @ cov @ grd)))
            rr = p / base_p
            glrr = g / p - base_g / base_p
            logrr_se = float(math.sqrt(max(0.0, glrr @ cov @ glrr)))
            local.append(
                {
                    "model": label, "interval_estimate": interval_label, "group": group,
                    "n": len(eval_data), "events": int(eval_data[outcome].sum()),
                    "adjusted_probability": p, "probability_lcl": max(0, p - 1.96 * se), "probability_ucl": min(1, p + 1.96 * se),
                    "risk_difference": rd, "rd_lcl": rd - 1.96 * rd_se, "rd_ucl": rd + 1.96 * rd_se,
                    "risk_ratio": rr, "rr_lcl": math.exp(math.log(rr) - 1.96 * logrr_se), "rr_ucl": math.exp(math.log(rr) + 1.96 * logrr_se),
                }
            )
        return local

    estimates = contrast_rows(z, "Pooled")
    if interaction:
        for intv in INTERVALS.values():
            eval_data = z.loc[z["interval"].eq(intv)].copy()
            eval_data["interval"] = intv
            estimates.extend(contrast_rows(eval_data, intv))
    coef = pd.DataFrame({"term": result.params.index, "estimate": result.params.values, "se": result.bse.values, "p": result.pvalues.values})
    coef["model"] = label
    interaction_terms = [i for i, name in enumerate(result.params.index) if ":" in name]
    interaction_p = np.nan
    if interaction_terms:
        beta = params[interaction_terms]
        vc = cov[np.ix_(interaction_terms, interaction_terms)]
        stat = float(beta @ np.linalg.pinv(vc) @ beta)
        interaction_p = float(chi2.sf(stat, len(interaction_terms)))
    return pd.DataFrame(estimates), coef, result, z, interaction_p


def observation_weights(base: pd.DataFrame) -> pd.DataFrame:
    z = base.copy()
    z["observed"] = z["observed_follow_profile"].astype(int)
    formula = (
        "observed ~ C(profile_base) + C(joint_group_base, Treatment(reference='Urban hukou/urban residence'))"
        " + C(interval, Treatment(reference='2011-2013')) + age_base + female_base"
        " + C(education_base, Treatment(reference='Lower')) + partnered_base"
        " + public_insurance_base + adl_difficulty_base + comorbidity_count_base"
    )
    cols = [
        "observed", "ID", "profile_base", "joint_group_base", "interval", "age_base",
        "female_base", "education_base", "partnered_base", "public_insurance_base",
        "adl_difficulty_base", "comorbidity_count_base", "weight_base",
    ]
    m = z[cols].dropna().copy()
    y, X = patsy.dmatrices(formula, m, return_type="dataframe")
    fit = sm.GEE(y.iloc[:, 0], X, groups=m["ID"].to_numpy(), family=sm.families.Binomial(), cov_struct=sm.cov_struct.Independence(), weights=m["weight_base"].to_numpy()).fit(maxiter=200)
    pred = np.clip(fit.predict(X), 0.05, 0.99)
    m["obs_prob"] = pred
    m["opw"] = 1.0 / pred
    lo, hi = m.loc[m.observed.eq(1), "opw"].quantile([0.01, 0.99])
    m["opw_trunc"] = m["opw"].clip(lo, hi)
    m["analysis_weight_opw"] = m["weight_base"] * m["opw_trunc"]
    return z.merge(m[["ID", "interval", "obs_prob", "opw_trunc", "analysis_weight_opw"]], on=["ID", "interval"], how="left")


def main():
    profiles = build_profiles()
    h = load_harmonized()
    pw = build_person_wave(h, profiles)
    intervals = build_intervals(pw)

    # Reproduce the original wave-specific samples before longitudinal modeling.
    wave_counts = []
    for w, year in WAVES.items():
        z = pw.loc[pw.wave.eq(w)]
        eligible = z.digest.eq(1) & z.profile.notna() & z.age.ge(45) & z.weight.gt(0)
        wave_counts.append({"wave": year, "eligible_n": int(eligible.sum()), "disease_positive_n": int((z.digest.eq(1) & z.age.ge(45)).sum())})
    pd.DataFrame(wave_counts).to_csv(OUT / "wave_sample_reproduction.csv", index=False)

    base = intervals.loc[intervals.eligible_base].copy()
    base["weight_norm"] = base.groupby("interval", observed=True).weight_base.transform(lambda x: x / x.mean())

    # Follow-up disposition by baseline profile and broad baseline state.
    disposition = np.select(
        [base.observed_follow_profile, base.known_death, base.alive_nonresponse, base.unknown_vital, base.responded_no_eligible_profile],
        ["Observed eligible profile", "Known death", "Alive nonresponse", "Vital status unknown", "Responded, no eligible profile"],
        default="Other/unclassified",
    )
    base["followup_disposition"] = disposition
    base["baseline_state2"] = np.where(base.none_listed_base, "None listed", "Any listed treatment")
    disp_rows = []
    for keys, g in base.groupby(["interval", "baseline_state2", "followup_disposition"], observed=True):
        denom = base.loc[(base.interval == keys[0]) & (base.baseline_state2 == keys[1]), "weight_norm"].sum()
        disp_rows.append({"interval": keys[0], "baseline_state": keys[1], "followup_disposition": keys[2], "n": len(g), "weighted_percent": 100 * g.weight_norm.sum() / denom})
    pd.DataFrame(disp_rows).to_csv(OUT / "followup_disposition_by_baseline_state.csv", index=False)

    # Diagnosis consistency and newly reported diagnosis.
    all_int = intervals.loc[intervals.age_base.ge(45) & intervals.weight_base.gt(0)].copy()
    diagnosis_rows = []
    for interval, g in all_int.groupby("interval", observed=True):
        both_resp = g.iwstat_base.eq(1) & g.iwstat_follow.eq(1)
        denom_yes = (both_resp & g.digest_base.eq(1)).sum()
        denom_no = (both_resp & g.digest_base.eq(0)).sum()
        diagnosis_rows.append(
            {
                "interval": interval,
                "both_interviewed_n": int(both_resp.sum()),
                "baseline_yes_n": int(denom_yes),
                "yes_to_no_n": int((both_resp & g.digest_base.eq(1) & g.digest_follow.eq(0)).sum()),
                "yes_to_no_percent": 100 * (both_resp & g.digest_base.eq(1) & g.digest_follow.eq(0)).sum() / denom_yes if denom_yes else np.nan,
                "baseline_no_n": int(denom_no),
                "no_to_yes_n": int((both_resp & g.digest_base.eq(0) & g.digest_follow.eq(1)).sum()),
                "no_to_yes_percent": 100 * (both_resp & g.digest_base.eq(0) & g.digest_follow.eq(1)).sum() / denom_no if denom_no else np.nan,
            }
        )
    pd.DataFrame(diagnosis_rows).to_csv(OUT / "diagnosis_consistency.csv", index=False)

    incident = all_int.loc[all_int.iwstat_base.eq(1) & all_int.iwstat_follow.eq(1) & all_int.digest_base.eq(0) & all_int.digest_follow.eq(1)].copy()
    incident_rows = []
    for (interval, profile), g in incident.groupby(["interval", "profile_follow"], dropna=False, observed=True):
        denom = incident.loc[incident.interval.eq(interval), "weight_base"].sum()
        incident_rows.append({"interval": interval, "followup_profile": profile if pd.notna(profile) else "Missing/ineligible", "n": len(g), "weighted_percent": 100 * g.weight_base.sum() / denom})
    pd.DataFrame(incident_rows).to_csv(OUT / "newly_reported_diagnosis_profiles.csv", index=False)

    observed = base.loc[base.observed_follow_profile].copy()
    transition4 = []
    for (interval, state), g in observed.groupby(["interval", "transition4"], observed=True):
        denom = observed.loc[observed.interval.eq(interval), "weight_norm"].sum()
        transition4.append({"interval": interval, "transition": state, "n": len(g), "weighted_percent": 100 * g.weight_norm.sum() / denom})
    pd.DataFrame(transition4).to_csv(OUT / "transition_four_type.csv", index=False)

    for interval, g in observed.groupby("interval", observed=True):
        weighted_table(g, "profile_base", "profile_follow", "weight_norm").to_csv(OUT / f"transition_matrix_{interval}.csv")

    # Counts used to describe the two baseline-state analysis cohorts.
    a = observed.loc[observed.none_listed_base].copy()
    a["initiation"] = a.any_listed_follow
    b = observed.loc[observed.any_listed_base].copy()
    b["cessation"] = b.none_listed_follow

    # Observation-probability weights for descriptive sensitivity estimates.
    opw = observation_weights(base)
    opw_obs = opw.loc[opw.observed_follow_profile & opw.analysis_weight_opw.notna()].copy()

    # Dispute-free sensitivity descriptive estimates.
    dispute_free = observed.loc[observed.dispute_base.eq(0) & observed.dispute_follow.eq(0)].copy()
    sensitivity_rows = []
    for source, d in [("Main", observed), ("Both dispute flags zero", dispute_free), ("Observation IPW", opw_obs)]:
        wt = "analysis_weight_opw" if source == "Observation IPW" else "weight_norm"
        for (interval, baseline_state), g in d.groupby(["interval", "baseline_state2"], observed=True):
            outcome = g.any_listed_follow if baseline_state == "None listed" else g.none_listed_follow
            sensitivity_rows.append({"analysis": source, "interval": interval, "baseline_state": baseline_state, "n": len(g), "weighted_event_percent": 100 * wmean(outcome.astype(float), g[wt])})
    pd.DataFrame(sensitivity_rows).to_csv(OUT / "transition_sensitivity_summary.csv", index=False)

    model_cols = [
        "ID", "joint_group_base", "interval", "age_base", "female_base",
        "education_base", "partnered_base", "public_insurance_base",
        "adl_difficulty_base", "comorbidity_count_base", "weight_norm",
    ]
    model_a = a[["initiation", *model_cols]].dropna()
    model_b = b[["cessation", *model_cols]].dropna()
    summary = {
        "wave_counts": wave_counts,
        "interval_baseline_n": {k: int(v) for k, v in base.groupby("interval", observed=True).size().items()},
        "observed_transition_n": {k: int(v) for k, v in observed.groupby("interval", observed=True).size().items()},
        "known_deaths": {k: int(v) for k, v in base.loc[base.known_death].groupby("interval", observed=True).size().items()},
        "eligible_outcome_a_n": int(len(a)),
        "eligible_outcome_b_n": int(len(b)),
        "complete_case_model_a_n": int(len(model_a)),
        "complete_case_model_a_events": int(model_a["initiation"].sum()),
        "complete_case_model_b_n": int(len(model_b)),
        "complete_case_model_b_events": int(model_b["cessation"].sum()),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
