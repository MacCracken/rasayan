//! AI-powered biochemistry queries via hoosh inference gateway.
//!
//! Provides a [`BiochemClient`] that wraps [`hoosh::HooshClient`] with
//! biochemistry-specific system prompts and query helpers.
//!
//! Requires the `ai` feature.
//!
//! # Example
//!
//! ```no_run
//! use rasayan::ai::BiochemClient;
//!
//! # async fn example() -> Result<(), Box<dyn std::error::Error>> {
//! let client = BiochemClient::new("http://localhost:8088", "llama3");
//! let answer = client.query("What is the Km of hexokinase for glucose?").await?;
//! println!("{}", answer.text);
//! # Ok(())
//! # }
//! ```

use hoosh::inference::{InferenceResponse, Message, Role};
use hoosh::{HooshClient, InferenceRequest};

/// System prompt establishing biochemistry expertise context.
const SYSTEM_PROMPT: &str = "\
You are a biochemistry expert assistant. Answer questions about enzyme kinetics, \
metabolic pathways, signal transduction, protein structure, membrane transport, \
and bioenergetics with precise, quantitative answers. Cite published Km, kcat, \
and Ki values when relevant. Use standard biochemistry nomenclature.";

/// A biochemistry-focused AI client wrapping hoosh.
#[derive(Clone)]
pub struct BiochemClient {
    client: HooshClient,
    model: String,
    system_prompt: String,
}

impl BiochemClient {
    /// Create a new client pointing at a hoosh instance.
    ///
    /// # Arguments
    /// - `base_url` — hoosh server URL (e.g., `"http://localhost:8088"`)
    /// - `model` — model identifier (e.g., `"llama3"`, `"claude-sonnet-4-20250514"`)
    #[must_use]
    pub fn new(base_url: impl Into<String>, model: impl Into<String>) -> Self {
        Self {
            client: HooshClient::new(base_url),
            model: model.into(),
            system_prompt: SYSTEM_PROMPT.to_string(),
        }
    }

    /// Override the default system prompt.
    #[must_use]
    pub fn with_system_prompt(mut self, prompt: impl Into<String>) -> Self {
        self.system_prompt = prompt.into();
        self
    }

    /// Send a single biochemistry query.
    ///
    /// Returns the full [`InferenceResponse`] including text, token usage,
    /// and latency.
    pub async fn query(&self, question: &str) -> Result<InferenceResponse, hoosh::HooshError> {
        let request = InferenceRequest {
            model: self.model.clone(),
            prompt: question.to_string(),
            system: Some(self.system_prompt.clone()),
            max_tokens: Some(2048),
            temperature: Some(0.3),
            ..Default::default()
        };
        self.client.infer(&request).await
    }

    /// Send a multi-turn conversation.
    ///
    /// Messages should alternate between [`Role::User`] and [`Role::Assistant`].
    /// The system prompt is automatically prepended.
    pub async fn chat(
        &self,
        messages: Vec<Message>,
    ) -> Result<InferenceResponse, hoosh::HooshError> {
        let mut all_messages = vec![Message::new(Role::System, &self.system_prompt)];
        all_messages.extend(messages);

        let request = InferenceRequest {
            model: self.model.clone(),
            messages: all_messages,
            max_tokens: Some(2048),
            temperature: Some(0.3),
            ..Default::default()
        };
        self.client.infer(&request).await
    }

    /// Ask about enzyme kinetics for a specific enzyme.
    pub async fn enzyme_query(
        &self,
        enzyme_name: &str,
    ) -> Result<InferenceResponse, hoosh::HooshError> {
        let prompt = format!(
            "Provide the published kinetic parameters for {enzyme_name}: \
             Km, kcat, catalytic efficiency (kcat/Km), known inhibitors, \
             optimal pH, and physiological function. Include units."
        );
        self.query(&prompt).await
    }

    /// Ask about a metabolic pathway.
    pub async fn pathway_query(
        &self,
        pathway_name: &str,
    ) -> Result<InferenceResponse, hoosh::HooshError> {
        let prompt = format!(
            "Describe the {pathway_name} pathway: key enzymes, substrates, products, \
             regulation points, ATP/NADH yield, and physiological significance. \
             Include rate-limiting steps."
        );
        self.query(&prompt).await
    }

    /// Ask about a protein sequence.
    pub async fn protein_query(
        &self,
        sequence: &str,
    ) -> Result<InferenceResponse, hoosh::HooshError> {
        let prompt = format!(
            "Analyze this protein sequence: {sequence}\n\
             Predict: likely function, domain composition, post-translational \
             modifications, subcellular localization, and any known homologs."
        );
        self.query(&prompt).await
    }

    /// Check if the hoosh server is reachable.
    pub async fn health(&self) -> Result<bool, hoosh::HooshError> {
        self.client.health().await
    }

    /// List available models on the hoosh server.
    pub async fn list_models(&self) -> Result<Vec<hoosh::inference::ModelInfo>, hoosh::HooshError> {
        self.client.list_models().await
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_client_creation() {
        let client = BiochemClient::new("http://localhost:8088", "llama3");
        assert_eq!(client.model, "llama3");
        assert!(client.system_prompt.contains("biochemistry"));
    }

    #[test]
    fn test_custom_system_prompt() {
        let client = BiochemClient::new("http://localhost:8088", "llama3")
            .with_system_prompt("Custom prompt");
        assert_eq!(client.system_prompt, "Custom prompt");
    }

    #[test]
    fn test_client_clone() {
        let client = BiochemClient::new("http://localhost:8088", "llama3");
        let _clone = client.clone();
    }
}
