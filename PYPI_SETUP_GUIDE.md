# PyPI Setup Guide for VCFX

This guide helps you set up PyPI publishing for the VCFX Python package.

## Method 1: Trusted Publishing (Recommended)

Trusted publishing uses OpenID Connect (OIDC) to securely publish without storing tokens.

### For TestPyPI

1. Go to https://test.pypi.org and log in
2. Navigate to your account settings
3. Click on "Publishing" in the sidebar
4. Add a new trusted publisher with these **exact** values:
   - **PyPI Project Name**: `vcfx` (or create the project first)
   - **Owner**: `jorgeMFS`
   - **Repository**: `VCFX`
   - **Workflow name**: `publish-pypi.yml`
   - **Environment**: `testpypi`

### For PyPI (Production)

1. Go to https://pypi.org and log in
2. Same process as above, but use:
   - **PyPI Project Name**: `vcfx`
   - **Owner**: `jorgeMFS`
   - **Repository**: `VCFX`
   - **Workflow name**: `publish-pypi.yml`
   - **Environment**: `pypi`

### Important Notes

- Repository name is **case-sensitive**
- The workflow file path must be `.github/workflows/publish-pypi.yml`
- The environment names must match exactly (`testpypi` and `pypi`)

## Method 2: API Token Authentication

If you need to publish immediately or prefer token-based auth:

### Step 1: Create API Tokens

**For TestPyPI:**
1. Go to https://test.pypi.org/manage/account/token/
2. Create a new API token
3. Set scope to "Entire account" or specific to `vcfx` project

**For PyPI:**
1. Go to https://pypi.org/manage/account/token/
2. Create a new API token
3. Set scope to "Entire account" or specific to `vcfx` project

### Step 2: Add Tokens to GitHub Secrets

1. Go to your repository settings on GitHub
2. Navigate to Secrets and variables → Actions
3. Add new repository secrets:
   - Name: `TEST_PYPI_API_TOKEN` → Value: Your TestPyPI token
   - Name: `PYPI_API_TOKEN` → Value: Your PyPI token

### Step 3: Update Workflow

Modify `.github/workflows/publish-pypi.yml`:

```yaml
# For TestPyPI - replace the publish step:
- name: Publish to TestPyPI
  uses: pypa/gh-action-pypi-publish@release/v1
  with:
    repository-url: https://test.pypi.org/legacy/
    user: __token__
    password: ${{ secrets.TEST_PYPI_API_TOKEN }}

# For PyPI - replace the publish step:
- name: Publish to PyPI
  uses: pypa/gh-action-pypi-publish@release/v1
  with:
    user: __token__
    password: ${{ secrets.PYPI_API_TOKEN }}
```

And remove the `permissions: id-token: write` section.

## Testing Your Setup

1. First test with TestPyPI:
   ```bash
   # Via GitHub Actions
   Go to Actions → Publish Python Package → Run workflow → Select "testpypi"
   
   # Or locally
   python scripts/build_and_publish.py --test
   ```

2. Verify installation:
   ```bash
   pip install -i https://test.pypi.org/simple/ vcfx
   ```

3. Once confirmed working, publish to PyPI:
   ```bash
   # Create a release on GitHub with tag v1.0.3
   # Or manually trigger the workflow selecting "pypi"
   ```

## Troubleshooting

### "invalid-publisher" Error
- Ensure the repository name matches exactly (case-sensitive)
- Check that the workflow filename is correct
- Verify the environment name matches

### "403 Forbidden" Error
- Token might be invalid or expired
- Token might not have permission for the project
- Project name might already be taken

### Package Name Already Exists
- If `vcfx` is taken, you'll need to choose a different name
- Update the package name in `python/pyproject.toml` 
