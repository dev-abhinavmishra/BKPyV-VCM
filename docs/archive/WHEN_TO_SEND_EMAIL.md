> **ARCHIVED - superseded historical record.** This document was generated during
> earlier development and contains stale or contradicted claims (test counts,
> parameter values, completion levels, citation details). Do not quote it.
> Canonical sources: ../README.md (root), ../BKPYV_MODEL_CARD.md,
> ../ISEF_PROJECT_OVERVIEW.md. Numbers must be regenerated from code, not
> copied from this file.

# Guide for Dr. Kotton Review - What to Send and When

## Files Generated for Clinical Review

1. **ONE_PAGE_SUMMARY_FOR_DR_KOTTON.md** - One-page overview (THIS IS THE MAIN DOCUMENT)
2. **outputs/clinical_review_plot.png** - Visual comparison of viral load trajectories
3. **outputs/clinical_review_summary.txt** - Detailed simulation metrics

## When to Send the Email

**Send the email when:**
- ✅ You have the three files listed above ready
- ✅ You have reviewed the plot and it looks reasonable to you
- ✅ You have tested the Streamlit dashboard and it runs without errors
- ✅ You feel confident explaining the project verbally

**DO NOT send yet if:**
- ❌ The simulation outputs look obviously wrong
- ❌ The dashboard has errors when you run it
- ❌ You don't understand what the plot is showing
- ❌ You can't explain the project basics in 2-3 minutes

## What to Send in the Email

**Email Body:** Use the email you drafted (it's excellent)

**Attachments:**
1. ONE_PAGE_SUMMARY_FOR_DR_KOTTON.md (rename to .pdf before sending if possible)
2. outputs/clinical_review_plot.png
3. outputs/clinical_review_summary.txt (optional)

**Optional Additional Attachments:**
- A screenshot of the Streamlit dashboard running
- A link to the GitHub repository (if it's public)

## What to Do Before Sending

### 1. Test the Dashboard
```bash
cd /Users/abhinavmishra/Projects/virtual-cell-model
python3 -m streamlit run src/vcm/ui/dashboard.py
```
- Does it load without errors?
- Can you see the different pages?
- Do the plots render correctly?

### 2. Review the Outputs
- Look at outputs/clinical_review_plot.png
- Does it make sense to you?
- Are the trends in the right direction (higher tacrolimus = higher viral load)?

### 3. Prepare for Questions
Dr. Kotton might ask:
- "How does the model work?" - Be ready to explain the pathway-based approach
- "What data is this based on?" - Mention the three validation sources
- "Have you validated this?" - Explain the external validation module
- "What makes this different from existing models?" - Explain the VCM approach

### 4. Convert One-Page Summary to PDF
```bash
# Option 1: Use Preview on Mac
# Open ONE_PAGE_SUMMARY_FOR_DR_KOTTON.md in Preview, File > Export as PDF

# Option 2: Use pandoc if installed
pandoc ONE_PAGE_SUMMARY_FOR_DR_KOTTON.md -o ONE_PAGE_SUMMARY_FOR_DR_KOTTON.pdf
```

## Email Sending Checklist

Before hitting send, verify:
- [ ] All attachments are included
- [ ] Your email address is in the contact section
- [ ] Your school name is correct
- [ ] You've proofread for typos
- [ ] You've run the dashboard successfully
- [ ] You can explain the project in under 5 minutes

## After Sending

**Wait for her response.** She may:
- Take a few days to review
- Forward to Dr. Hirsch (since she already referred you)
- Ask for more details
- Offer to meet or have a call

**If she agrees to the review:**
- Be prepared to send more detailed simulation results
- Consider scheduling a brief call if she offers
- Be ready to make adjustments based on her feedback

**If she asks for a demo:**
- Have the Streamlit dashboard ready to run
- Practice showing her the key features:
  - Risk prediction page
  - Drug tapering simulation
  - Sensitivity analysis

## Backup Plan

If she doesn't respond within 2 weeks:
- Send a polite follow-up email
- Focus on completing other aspects of the project
- Consider reaching out to other mentors (Dr. Hirsch, local professors)

Remember: **Clinical plausibility review is a bonus, not a requirement** for ISEF. The project already has solid literature-based validation.
