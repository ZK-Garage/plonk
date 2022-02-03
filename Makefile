help: ## Display this help screen
	@grep -h -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'

doc: ## Generate documentation
	@cargo rustdoc --lib -- -D warnings

doc-internal: ## Generate documentation with private items
	@cargo rustdoc --lib -- --document-private-items -D warnings

book: ## Generate book
	@mdbook build plonk-book
	

.PHONY: help doc doc-internal book
